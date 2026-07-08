function prevPow2(n){var p=1;while((p<<1)<=n)p<<=1;return p}
function analyzeSpectrumFixed(){
  if(!st){message('Build atom first, then run the simulation before spectrum analysis.');return}
  var count=st.records?st.records.length:0;
  if(count<32){message('not enough probe records for FFT: '+count+'/32. Press Start or Atom oscillation only and let it run first.');return}
  var recs=st.records;
  var n=Math.min(prevPow2(recs.length),8192);
  if(n<32){message('not enough records after FFT window selection');return}
  var start=recs.length-n;
  var dt=recs[start+1].t-recs[start].t;
  if(!(dt>0)){message('invalid probe time step; clear records and run again');return}
  var v1=st.md.vLight[0],v2=st.md.vLight[1];
  var left=new Float64Array(n),right=new Float64Array(n),l2=new Float64Array(n),r2=new Float64Array(n);
  for(var i=0;i<n;i++){
    var r=recs[start+i];
    var win=0.5*(1-Math.cos(2*Math.PI*i/(n-1)));
    left[i]=win*(v1*r.dphi1_L1+v2*r.dphi2_L1);
    right[i]=win*(v1*r.dphi1_R1+v2*r.dphi2_R1);
    l2[i]=win*(v1*r.dphi1_L2+v2*r.dphi2_L2);
    r2[i]=win*(v1*r.dphi1_R2+v2*r.dphi2_R2);
  }
  var imL=new Float64Array(n),imR=new Float64Array(n),imL2=new Float64Array(n),imR2=new Float64Array(n);
  fftComplex(left,imL);fftComplex(right,imR);fftComplex(l2,imL2);fftComplex(r2,imR2);
  var half=n/2,omega=[],leftPower=[],rightPower=[],leftMovingLeft=[],rightMovingLeft=[],leftMovingRight=[],rightMovingRight=[];
  var dxProbe=st.p.probeSep;
  for(var k=1;k<half;k++){
    var om=2*Math.PI*k/(n*dt);omega.push(om);
    var pL=left[k]*left[k]+imL[k]*imL[k],pR=right[k]*right[k]+imR[k]*imR[k];
    leftPower.push(pL);rightPower.push(pR);
    var kk=Math.sqrt(Math.max(om*om-st.md.mLight2,0));
    if(kk*dxProbe<1e-8){leftMovingLeft.push(0);rightMovingLeft.push(pL);leftMovingRight.push(0);rightMovingRight.push(pR);continue}
    var c=Math.cos(kk*dxProbe),s=Math.sin(kk*dxProbe),det=2*s;
    if(Math.abs(det)<1e-8){leftMovingLeft.push(0);rightMovingLeft.push(pL);leftMovingRight.push(0);rightMovingRight.push(pR);continue}
    function split(Are,Aim,Bre,Bim){
      var Rre=(Bim-Aim*c-Are*s)/det;
      var Rim=(-Bre+Are*c-Aim*s)/det;
      var Lre=Are-Rre,Lim=Aim-Rim;
      return[Lre*Lre+Lim*Lim,Rre*Rre+Rim*Rim]
    }
    var spL=split(left[k],imL[k],l2[k],imL2[k]);
    var spR=split(right[k],imR[k],r2[k],imR2[k]);
    leftMovingLeft.push(spL[0]);rightMovingLeft.push(spL[1]);leftMovingRight.push(spR[0]);rightMovingRight.push(spR[1]);
  }
  st.lastSpectrum={omega:omega,leftPower:leftPower,rightPower:rightPower,leftMovingLeft:leftMovingLeft,rightMovingLeft:rightMovingLeft,leftMovingRight:leftMovingRight,rightMovingRight:rightMovingRight};
  drawSpectrum(st.lastSpectrum);drawDirection(st.lastSpectrum);
  var box=document.getElementById('peakReadout');
  if(box&&typeof spectrumPeakReport==='function')box.textContent=spectrumPeakReport(st.lastSpectrum);
  message('spectrum computed from '+n+' records, t window '+(recs[start].t).toFixed(2)+' to '+(recs[recs.length-1].t).toFixed(2));
}
(function(){
  var btn=document.getElementById('spectrumBtn');
  if(btn){btn.addEventListener('click',function(e){e.preventDefault();e.stopImmediatePropagation();analyzeSpectrumFixed();},true)}
  var peak=document.getElementById('peakBtn');
  if(peak){peak.addEventListener('click',function(e){e.preventDefault();e.stopImmediatePropagation();var box=document.getElementById('peakReadout');if(box)box.textContent=(typeof spectrumPeakReport==='function')?spectrumPeakReport(st&&st.lastSpectrum):'peak finder unavailable';},true)}
})();
