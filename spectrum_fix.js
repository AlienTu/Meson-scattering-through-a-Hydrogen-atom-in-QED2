function safeMsg(txt){if(typeof message==='function')message(txt);else console.log(txt)}
function simpleSpectrumDFT(){
  try{
    if(!st){safeMsg('Build atom first.');return}
    var recs=st.records||[];
    var count=recs.length;
    if(count<8){safeMsg('not enough probe records for Fourier transform: '+count+'/8. Press Start or Atom oscillation only and let the trace run.');return}
    var n=Math.min(count,4096);
    var start=count-n;
    var t0=recs[start].t;
    var t1=recs[count-1].t;
    var T=t1-t0;
    if(!(T>0)){safeMsg('invalid time window for spectrum; clear records and run again.');return}
    var v1=st.md.vLight[0],v2=st.md.vLight[1];
    var yL=new Float64Array(n),yR=new Float64Array(n),yL2=new Float64Array(n),yR2=new Float64Array(n);
    var meanL=0,meanR=0,meanL2=0,meanR2=0;
    for(var i=0;i<n;i++){
      var r=recs[start+i];
      var aL=v1*r.dphi1_L1+v2*r.dphi2_L1;
      var aR=v1*r.dphi1_R1+v2*r.dphi2_R1;
      var aL2=v1*r.dphi1_L2+v2*r.dphi2_L2;
      var aR2=v1*r.dphi1_R2+v2*r.dphi2_R2;
      yL[i]=aL;yR[i]=aR;yL2[i]=aL2;yR2[i]=aR2;
      meanL+=aL;meanR+=aR;meanL2+=aL2;meanR2+=aR2;
    }
    meanL/=n;meanR/=n;meanL2/=n;meanR2/=n;
    var omegaMax=12;
    var bins=520;
    var omega=[],leftPower=[],rightPower=[],leftMovingLeft=[],rightMovingLeft=[],leftMovingRight=[],rightMovingRight=[];
    var dxProbe=st.p.probeSep;
    function coeff(y,mean,om){
      var re=0,im=0,wsum=0;
      for(var i=0;i<n;i++){
        var tau=recs[start+i].t-t0;
        var win=0.5*(1-Math.cos(2*Math.PI*i/Math.max(n-1,1)));
        var yy=(y[i]-mean)*win;
        var ph=om*tau;
        re+=yy*Math.cos(ph);
        im-=yy*Math.sin(ph);
        wsum+=win;
      }
      var sc=wsum>0?1/wsum:1;
      return[re*sc,im*sc];
    }
    for(var b=1;b<=bins;b++){
      var om=omegaMax*b/bins;
      omega.push(om);
      var L=coeff(yL,meanL,om),R=coeff(yR,meanR,om),L2=coeff(yL2,meanL2,om),R2=coeff(yR2,meanR2,om);
      var pL=L[0]*L[0]+L[1]*L[1],pR=R[0]*R[0]+R[1]*R[1];
      leftPower.push(pL);rightPower.push(pR);
      var kk=Math.sqrt(Math.max(om*om-st.md.mLight2,0));
      if(kk*dxProbe<1e-8){leftMovingLeft.push(0);rightMovingLeft.push(pL);leftMovingRight.push(0);rightMovingRight.push(pR);continue}
      var c=Math.cos(kk*dxProbe),s=Math.sin(kk*dxProbe),det=2*s;
      if(Math.abs(det)<1e-8){leftMovingLeft.push(0);rightMovingLeft.push(pL);leftMovingRight.push(0);rightMovingRight.push(pR);continue}
      function split(A,B){
        var Are=A[0],Aim=A[1],Bre=B[0],Bim=B[1];
        var Rre=(Bim-Aim*c-Are*s)/det;
        var Rim=(-Bre+Are*c-Aim*s)/det;
        var Lre=Are-Rre,Lim=Aim-Rim;
        return[Lre*Lre+Lim*Lim,Rre*Rre+Rim*Rim]
      }
      var spL=split(L,L2),spR=split(R,R2);
      leftMovingLeft.push(spL[0]);rightMovingLeft.push(spL[1]);leftMovingRight.push(spR[0]);rightMovingRight.push(spR[1]);
    }
    st.lastSpectrum={omega:omega,leftPower:leftPower,rightPower:rightPower,leftMovingLeft:leftMovingLeft,rightMovingLeft:rightMovingLeft,leftMovingRight:leftMovingRight,rightMovingRight:rightMovingRight};
    if(typeof drawSpectrum==='function')drawSpectrum(st.lastSpectrum);
    if(typeof drawDirection==='function')drawDirection(st.lastSpectrum);
    var box=document.getElementById('peakReadout');
    if(box&&typeof spectrumPeakReport==='function')box.textContent=spectrumPeakReport(st.lastSpectrum);
    safeMsg('direct Fourier spectrum computed from '+n+' records, time window '+t0.toFixed(2)+' to '+t1.toFixed(2));
  }catch(err){
    safeMsg('spectrum error: '+(err&&err.message?err.message:String(err)));
    console.error(err);
  }
}
(function(){
  function bindSpectrum(){
    var btn=document.getElementById('spectrumBtn');
    if(btn){btn.onclick=null;btn.addEventListener('click',function(e){e.preventDefault();e.stopImmediatePropagation();simpleSpectrumDFT();},true)}
    var peak=document.getElementById('peakBtn');
    if(peak){peak.onclick=null;peak.addEventListener('click',function(e){e.preventDefault();e.stopImmediatePropagation();var box=document.getElementById('peakReadout');if(box)box.textContent=(typeof spectrumPeakReport==='function')?spectrumPeakReport(st&&st.lastSpectrum):'peak finder unavailable';},true)}
  }
  if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',bindSpectrum);else bindSpectrum();
})();
