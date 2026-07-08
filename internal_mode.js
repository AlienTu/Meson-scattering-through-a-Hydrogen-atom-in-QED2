function imDot(a,b){var s=0;for(var i=0;i<a.length;i++)s+=a[i]*b[i];return s}
function imNorm(a){return Math.sqrt(Math.max(imDot(a,a),0))}
function imNormalize(a){var n=imNorm(a);if(n>0)for(var i=0;i<a.length;i++)a[i]/=n;return a}
function imAxpy(y,a,x){for(var i=0;i<y.length;i++)y[i]+=a*x[i]}
function imScale(y,a){for(var i=0;i<y.length;i++)y[i]*=a}

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
function jacobiSmall(A,maxSweeps,tol){
  var n=A.length,V=[];for(var i=0;i<n;i++){V[i]=new Float64Array(n);V[i][i]=1}
  for(var sw=0;sw<maxSweeps;sw++){
    var p=0,q=1,mx=0;for(var a=0;a<n;a++)for(var b=a+1;b<n;b++){var val=Math.abs(A[a][b]);if(val>mx){mx=val;p=a;q=b}}
    if(mx<tol)break;var app=A[p][p],aqq=A[q][q],apq=A[p][q];if(Math.abs(apq)<1e-300)continue;
    var tau=(aqq-app)/(2*apq),t=Math.sign(tau||1)/(Math.abs(tau)+Math.sqrt(1+tau*tau)),c=1/Math.sqrt(1+t*t),s=t*c;
    for(var k=0;k<n;k++)if(k!==p&&k!==q){var akp=A[k][p],akq=A[k][q];A[k][p]=A[p][k]=c*akp-s*akq;A[k][q]=A[q][k]=s*akp+c*akq}
    A[p][p]=c*c*app-2*s*c*apq+s*s*aqq;A[q][q]=s*s*app+2*s*c*apq+c*c*aqq;A[p][q]=A[q][p]=0;
    for(var kk=0;kk<n;kk++){var vkp=V[kk][p],vkq=V[kk][q];V[kk][p]=c*vkp-s*vkq;V[kk][q]=s*vkp+c*vkq}
  }
  var pairs=[];for(var ii=0;ii<n;ii++)pairs.push({value:A[ii][ii],index:ii});pairs.sort(function(a,b){return a.value-b.value});return{pairs:pairs,V:V};
}
function lanczosModes(op,kdim){
  var dim=2*op.n,k=Math.min(kdim,dim-2),qPrev=new Float64Array(dim),q=new Float64Array(dim),z=new Float64Array(dim),Q=[],alpha=[],beta=[],b=0;
  for(var i=0;i<dim;i++){var site=Math.floor(i/2),x=op.x[site]||0,env=Math.exp(-0.5*x*x/225);q[i]=env*(Math.sin(0.173*(i+1))+0.37*Math.cos(0.117*(i+1)))}
  imNormalize(q);
  for(var j=0;j<k;j++){
    Q.push(new Float64Array(q));z.fill(0);applyFluctuationOp(q,op,z);if(j>0)imAxpy(z,-b,qPrev);var a=imDot(q,z);alpha.push(a);imAxpy(z,-a,q);
    for(var r=0;r<Q.length;r++){var c=imDot(Q[r],z);imAxpy(z,-c,Q[r])}
    b=imNorm(z);if(j<k-1)beta.push(b);if(b<1e-12)break;qPrev.set(q);q.set(z);imScale(q,1/b);
  }
  var m=alpha.length,T=[];for(var t=0;t<m;t++){T[t]=new Float64Array(m);T[t][t]=alpha[t];if(t<m-1){T[t][t+1]=beta[t];T[t+1][t]=beta[t]}}
  return{Q:Q,eig:jacobiSmall(T,160,1e-10)};
}
function reconstructLanczosVector(Q,V,col){
  var dim=Q[0].length,u=new Float64Array(dim);
  for(var j=0;j<Q.length;j++){var c=V[j][col];for(var i=0;i<dim;i++)u[i]+=c*Q[j][i]}
  imNormalize(u);return u;
}
function analyzeModeFromVector(v,op,lambda){
  var n=op.n,dx=op.dx,rho=new Float64Array(n),nm=0,hv=new Float64Array(2*n),tmp=new Float64Array(2*n);
  applyFluctuationOp(v,op,hv);var ray=imDot(v,hv);for(var i=0;i<2*n;i++)tmp[i]=hv[i]-ray*v[i];
  for(var r=0;r<n;r++){var den=v[2*r]*v[2*r]+v[2*r+1]*v[2*r+1];rho[r]=den;nm+=den*dx}
  for(var q=0;q<n;q++)rho[q]/=Math.max(nm,1e-300);
  var ipr=0,wR=0,edge=0,R=Math.min(18,0.25*st.p.L),edgeR=0.42*st.p.L;
  for(var k=0;k<n;k++){ipr+=rho[k]*rho[k]*dx;if(Math.abs(op.x[k])<R)wR+=rho[k]*dx;if(Math.abs(op.x[k])>edgeR)edge+=rho[k]*dx}
  return{lambda:ray,omega:Math.sqrt(Math.max(ray,0)),vec:v,rho:rho,ipr:ipr,wR:wR,edge:edge,relres:imNorm(tmp)/Math.max(1,Math.abs(ray))};
}
function fullVectorFromInterleaved(v,op){
  var N=st.p.N,out=new Float64Array(2*N);
  for(var i=0;i<op.n;i++){var site=i+1;out[site]=v[2*i];out[N+site]=v[2*i+1]}
  var dx=st.p.dx,s=0;for(var k=0;k<out.length;k++)s+=out[k]*out[k]*dx;s=Math.sqrt(Math.max(s,1e-300));for(var q=0;q<out.length;q++)out[q]/=s;out[0]=0;out[N-1]=0;out[N]=0;out[2*N-1]=0;return out;
}
function solveInternalMode(){
  if(!st)buildAtom();var op=buildFluctuationOp(),lo=st.md.mLight,kdim=90,lan=lanczosModes(op,kdim),modes=[];
  for(var p=0;p<Math.min(40,lan.eig.pairs.length);p++){
    var pair=lan.eig.pairs[p],v=reconstructLanczosVector(lan.Q,lan.eig.V,pair.index),m=analyzeModeFromVector(v,op,pair.value);
    if(m.omega<1.15*lo&&m.relres<5e-4)modes.push(m);
  }
  modes.sort(function(a,b){return a.lambda-b.lambda});
  var chosen=null,zeroCut=0.25;
  for(var i=0;i<modes.length;i++){var x=modes[i];if(x.omega>zeroCut&&x.omega<lo&&x.wR>0.45&&x.edge<0.05){chosen=x;break}}
  if(!chosen){for(var j=0;j<modes.length;j++){var y=modes[j];if(y.omega>zeroCut&&y.omega<lo&&y.wR>0.35){chosen=y;break}}}
  if(!chosen&&modes.length)chosen=modes[0];
  if(!chosen){var b0=document.getElementById('internalModeReadout');if(b0)b0.textContent='No acceptable mode found. Check static background.';return null}
  st.internalMode={u:fullVectorFromInterleaved(chosen.vec,op),omega:chosen.omega,lambda:chosen.lambda,candidates:modes};
  var lines=[];lines.push('full-space Lanczos eigensolver; no parity projection');lines.push('static background still uses half-axis symmetry, but modes are solved on full x-space');lines.push('threshold m_light = '+lo.toFixed(6));lines.push('chosen mode: omega = '+chosen.omega.toFixed(6)+', omega^2 = '+chosen.lambda.toFixed(6));lines.push('center weight = '+chosen.wR.toFixed(4)+', edge weight = '+chosen.edge.toExponential(2)+', IPR = '+chosen.ipr.toFixed(4));lines.push('relative residual = '+chosen.relres.toExponential(3));lines.push('low candidate modes:');
  for(var r=0;r<Math.min(14,modes.length);r++){var cc=modes[r];lines.push(r+': omega='+cc.omega.toFixed(6)+', wR='+cc.wR.toFixed(3)+', edge='+cc.edge.toExponential(1)+', res='+cc.relres.toExponential(1))}
  var box=document.getElementById('internalModeReadout');if(box)box.textContent=lines.join('\n');if(typeof message==='function')message('internal mode solved on full space: omega = '+chosen.omega.toFixed(4));return st.internalMode;
}
function applyInternalMode(){
  if(!st||!st.internalMode)solveInternalMode();if(!st||!st.internalMode)return;var ampEl=document.getElementById('internalAmp'),phaseEl=document.getElementById('internalPhase'),amp=ampEl?parseFloat(ampEl.value):0.08,phase=phaseEl?parseFloat(phaseEl.value):0,p=st.p,n=p.N,u=st.internalMode.u,omega=st.internalMode.omega;
  st.running=false;st.phi1=st.phi1Static.slice();st.phi2=st.phi2Static.slice();st.pi1=new Float64Array(n);st.pi2=new Float64Array(n);st.records=[];st.lastSpectrum=null;st.time=0;st.step=0;
  for(var i=1;i<n-1;i++){st.phi1[i]+=amp*Math.cos(phase)*u[i];st.phi2[i]+=amp*Math.cos(phase)*u[n+i];st.pi1[i]+=-amp*omega*Math.sin(phase)*u[i];st.pi2[i]+=-amp*omega*Math.sin(phase)*u[n+i]}
  st.phi1[0]=st.phi1Static[0];st.phi2[0]=st.phi2Static[0];st.phi1[n-1]=st.phi1Static[n-1];st.phi2[n-1]=st.phi2Static[n-1];
  if(p.driveMode==='packet')addIncomingPacket(st);st.initialEnergy=totalEnergy(st);drawAll();updateReadouts();if(typeof message==='function')message('internal vibration applied');
}
function atomOnlyRun(){var dm=document.getElementById('driveMode');if(dm)dm.value='none';applyInternalMode();if(st)st.running=true}
function setupInternalButtons(){var solveBtn=document.getElementById('solveInternalBtn'),applyBtn=document.getElementById('applyInternalBtn'),atomBtn=document.getElementById('atomOnlyBtn');if(solveBtn)solveBtn.addEventListener('click',solveInternalMode);if(applyBtn)applyBtn.addEventListener('click',applyInternalMode);if(atomBtn)atomBtn.addEventListener('click',atomOnlyRun)}
setupInternalButtons();
