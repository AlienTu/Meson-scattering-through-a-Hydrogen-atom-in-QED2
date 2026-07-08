function imDot(a,b){var s=0;for(var i=0;i<a.length;i++)s+=a[i]*b[i];return s}
function imNorm(a){return Math.sqrt(Math.max(imDot(a,a),0))}
function imNormalize(a){var n=imNorm(a);if(n>0)for(var i=0;i<a.length;i++)a[i]/=n;return a}
function imAxpy(y,a,x){for(var i=0;i<y.length;i++)y[i]+=a*x[i]}

function buildFluctuationOp(){
  var p=st.p,n=p.N-2,dx=p.dx,inv=1/(dx*dx),beta=TWO_SQRT_PI,mix=p.g*p.g;
  var D0=new Float64Array(n),D1=new Float64Array(n),D2=new Float64Array(n),xs=new Float64Array(n);
  for(var i=0;i<n;i++){
    var j=i+1;xs[i]=st.x[j];
    D0[i]=2*inv+mix+2*Math.PI*p.k1*Math.cos(beta*st.phi1Static[j]);
    D1[i]=mix;
    D2[i]=2*inv+mix+2*Math.PI*p.k2*Math.cos(beta*st.phi2Static[j]);
  }
  return{n:n,dx:dx,inv:inv,D0:D0,D1:D1,D2:D2,x:xs};
}
function applyFluctuationOp(v,op,out){
  var n=op.n,inv=op.inv;
  for(var i=0;i<n;i++){
    var a=2*i,b=a+1,u1=v[a],u2=v[b];
    var l1=i>0?v[a-2]:0,r1=i<n-1?v[a+2]:0,l2=i>0?v[a-1]:0,r2=i<n-1?v[a+3]:0;
    out[a]=op.D0[i]*u1+op.D1[i]*u2-inv*(l1+r1);
    out[b]=op.D1[i]*u1+op.D2[i]*u2-inv*(l2+r2);
  }
}
function inv2sym(a,b,c){var det=a*c-b*b;if(Math.abs(det)<1e-13)det=det>=0?1e-13:-1e-13;return[c/det,-b/det,a/det]}
function solveShiftedOp(op,sigma,b){
  var n=op.n,inv=op.inv,coef=inv*inv,E0=new Float64Array(n),E1=new Float64Array(n),E2=new Float64Array(n),G0=new Float64Array(n),G1=new Float64Array(n),prev=null;
  for(var i=0;i<n;i++){
    var a=op.D0[i]-sigma,c=op.D2[i]-sigma,bb=op.D1[i],r0=b[2*i],r1=b[2*i+1];
    if(i>0){a-=coef*prev[0];bb-=coef*prev[1];c-=coef*prev[2];r0+=inv*G0[i-1];r1+=inv*G1[i-1]}
    var iv=inv2sym(a,bb,c);G0[i]=iv[0]*r0+iv[1]*r1;G1[i]=iv[1]*r0+iv[2]*r1;E0[i]=-inv*iv[0];E1[i]=-inv*iv[1];E2[i]=-inv*iv[2];prev=iv;
  }
  var x=new Float64Array(2*n);
  for(var j=n-1;j>=0;j--){var x0=G0[j],x1=G1[j];if(j<n-1){var y0=x[2*j+2],y1=x[2*j+3];x0-=E0[j]*y0+E1[j]*y1;x1-=E1[j]*y0+E2[j]*y1}x[2*j]=x0;x[2*j+1]=x1}
  return x;
}
function projectParity(v,op,par){
  var n=op.n;
  for(var i=0;i<Math.floor(n/2);i++){
    var j=n-1-i,a=2*i,b=a+1,c=2*j,d=c+1,p1=0.5*(v[a]+par*v[c]),p2=0.5*(v[b]+par*v[d]);
    v[a]=p1;v[b]=p2;v[c]=par*p1;v[d]=par*p2;
  }
  if(n%2===1&&par<0){var m=Math.floor(n/2);v[2*m]=0;v[2*m+1]=0}
  return v;
}
function parityError(v,op,par){
  var n=op.n,e=0,s=0;
  for(var i=0;i<n;i++){var j=n-1-i,a=2*i,b=a+1,c=2*j,d=c+1,d1=v[a]-par*v[c],d2=v[b]-par*v[d];e+=d1*d1+d2*d2;s+=v[a]*v[a]+v[b]*v[b]}
  return Math.sqrt(e/Math.max(s,1e-300));
}
function orthogonalize(v,vecs){for(var k=0;k<vecs.length;k++){var u=vecs[k],c=imDot(v,u);for(var i=0;i<v.length;i++)v[i]-=c*u[i]}}
function inverseMode(op,sigma,prev,iters,par){
  var dim=2*op.n,v=new Float64Array(dim),hv=new Float64Array(dim),tmp=new Float64Array(dim);
  for(var i=0;i<dim;i++)v[i]=Math.sin(0.37*(i+1)+sigma)+0.31*Math.cos(0.119*(i+2)+0.7*sigma);
  projectParity(v,op,par);orthogonalize(v,prev);imNormalize(v);
  for(var it=0;it<iters;it++){var y=solveShiftedOp(op,sigma,v);projectParity(y,op,par);orthogonalize(y,prev);projectParity(y,op,par);imNormalize(y);v=y}
  applyFluctuationOp(v,op,hv);var lam=imDot(v,hv);for(var j=0;j<dim;j++)tmp[j]=hv[j]-lam*v[j];
  return{lambda:lam,omega:Math.sqrt(Math.max(lam,0)),vec:v,res:imNorm(tmp),relres:imNorm(tmp)/Math.max(1,Math.abs(lam)),parity:par>0?'even':'odd',perr:parityError(v,op,par)};
}
function analyzeMode(m,op,R){
  var n=op.n,dx=op.dx,rho=new Float64Array(n),nm=0;
  for(var i=0;i<n;i++){var r=m.vec[2*i]*m.vec[2*i]+m.vec[2*i+1]*m.vec[2*i+1];rho[i]=r;nm+=r*dx}
  for(var j=0;j<n;j++)rho[j]/=Math.max(nm,1e-300);
  var ipr=0,wR=0;for(var k=0;k<n;k++){ipr+=rho[k]*rho[k]*dx;if(Math.abs(op.x[k])<R)wR+=rho[k]*dx}
  m.rho=rho;m.ipr=ipr;m.wR=wR;return m;
}
function mergeModes(ms){
  ms.sort(function(a,b){return a.lambda-b.lambda});var out=[];
  for(var i=0;i<ms.length;i++){var m=ms[i];if(!Number.isFinite(m.lambda))continue;var last=out[out.length-1];if(last&&last.parity===m.parity&&Math.abs(m.lambda-last.lambda)<2e-4){if(m.relres<last.relres)out[out.length-1]=m}else out.push(m)}
  return out;
}
function fullVectorFromInterleaved(v,op){
  var N=st.p.N,out=new Float64Array(2*N);
  for(var i=0;i<op.n;i++){var site=i+1;out[site]=v[2*i];out[N+site]=v[2*i+1]}
  var dx=st.p.dx,s=0;for(var k=0;k<out.length;k++)s+=out[k]*out[k]*dx;s=Math.sqrt(Math.max(s,1e-300));for(var q=0;q<out.length;q++)out[q]/=s;return out;
}
function solveInternalMode(){
  if(!st)buildAtom();var op=buildFluctuationOp(),lo=st.md.mLight,lo2=lo*lo,R=Math.min(18,0.25*st.p.L),iters=10;
  var shifts=[];shifts.push(-0.02,0.02,0.05*lo2,0.10*lo2,0.20*lo2,0.35*lo2,0.50*lo2,0.65*lo2,0.80*lo2,0.92*lo2,0.98*lo2);
  for(var j=0;j<12;j++)shifts.push(lo2*(j+1)/13);
  var modes=[];
  for(var pidx=0;pidx<2;pidx++){
    var par=pidx===0?1:-1,prev=[];
    for(var s=0;s<shifts.length;s++){
      var m=inverseMode(op,shifts[s],prev,iters,par);m=analyzeMode(m,op,R);
      if(m.omega<1.15*lo&&m.relres<1e-5){modes.push(m);prev.push(m.vec)}
      if(prev.length>6)prev.shift();
    }
  }
  modes=mergeModes(modes);
  modes.sort(function(a,b){return a.lambda-b.lambda});
  var chosen=null,zeroCut=0.25;
  for(var i=0;i<modes.length;i++){var x=modes[i];if(x.omega>zeroCut&&x.omega<lo&&x.wR>0.45){chosen=x;break}}
  if(!chosen){for(var ii=0;ii<modes.length;ii++){var y=modes[ii];if(y.omega>zeroCut&&y.wR>0.35){chosen=y;break}}}
  if(!chosen&&modes.length)chosen=modes[0];
  if(!chosen){var b0=document.getElementById('internalModeReadout');if(b0)b0.textContent='No mode found. Check static background.';return null}
  st.internalMode={u:fullVectorFromInterleaved(chosen.vec,op),omega:chosen.omega,lambda:chosen.lambda,candidates:modes};
  var lines=[];lines.push('method copied from Hydrogen-formation-quench spectrum solver');lines.push('block-tridiagonal operator in interleaved (phi1, phi2) basis');lines.push('gauge convention: M12 = g^2, threshold m_light = '+lo.toFixed(6));lines.push('chosen mode: omega = '+chosen.omega.toFixed(6)+', omega^2 = '+chosen.lambda.toFixed(6));lines.push('parity = '+chosen.parity+', center weight = '+chosen.wR.toFixed(4)+', IPR = '+chosen.ipr.toFixed(4));lines.push('relative residual = '+chosen.relres.toExponential(3)+', parity error = '+chosen.perr.toExponential(2));lines.push('low candidate modes:');
  for(var r=0;r<Math.min(12,modes.length);r++){var cc=modes[r];lines.push(r+': omega='+cc.omega.toFixed(6)+', parity='+cc.parity+', wR='+cc.wR.toFixed(3)+', res='+cc.relres.toExponential(1))}
  var box=document.getElementById('internalModeReadout');if(box)box.textContent=lines.join('\n');if(typeof message==='function')message('internal mode solved by reference-style block solver: omega = '+chosen.omega.toFixed(4));return st.internalMode;
}
function applyInternalMode(){
  if(!st||!st.internalMode)solveInternalMode();if(!st||!st.internalMode)return;var ampEl=document.getElementById('internalAmp'),phaseEl=document.getElementById('internalPhase'),amp=ampEl?parseFloat(ampEl.value):0.08,phase=phaseEl?parseFloat(phaseEl.value):0,p=st.p,n=p.N,u=st.internalMode.u,omega=st.internalMode.omega;
  st.running=false;st.phi1=st.phi1Static.slice();st.phi2=st.phi2Static.slice();st.pi1=new Float64Array(n);st.pi2=new Float64Array(n);st.records=[];st.lastSpectrum=null;st.time=0;st.step=0;
  for(var i=1;i<n-1;i++){st.phi1[i]+=amp*Math.cos(phase)*u[i];st.phi2[i]+=amp*Math.cos(phase)*u[n+i];st.pi1[i]+=-amp*omega*Math.sin(phase)*u[i];st.pi2[i]+=-amp*omega*Math.sin(phase)*u[n+i]}
  if(p.driveMode==='packet')addIncomingPacket(st);st.initialEnergy=totalEnergy(st);drawAll();updateReadouts();if(typeof message==='function')message('internal vibration applied');
}
function atomOnlyRun(){var dm=document.getElementById('driveMode');if(dm)dm.value='none';applyInternalMode();if(st)st.running=true}
function setupInternalButtons(){var solveBtn=document.getElementById('solveInternalBtn'),applyBtn=document.getElementById('applyInternalBtn'),atomBtn=document.getElementById('atomOnlyBtn');if(solveBtn)solveBtn.addEventListener('click',solveInternalMode);if(applyBtn)applyBtn.addEventListener('click',applyInternalMode);if(atomBtn)atomBtn.addEventListener('click',atomOnlyRun)}
setupInternalButtons();
