function hqDot(a,b){var s=0;for(var i=0;i<a.length;i++)s+=a[i]*b[i];return s}
function hqNorm(a){return Math.sqrt(Math.max(hqDot(a,a),0))}
function hqApplyJHalf(v,out,p1,p2,g,k1,k2,dx){
  var sp=Math.sqrt(Math.PI), beta=2*sp, Nh=p1.length, n=Nh-2, inv=1/(dx*dx), mix=g*g;
  for(var i=0;i<n;i++){
    var gi=i+1,u1=v[i],u2=v[i+n];
    var l1=i>0?v[i-1]:0,r1=i<n-1?v[i+1]:0,l2=i>0?v[i+n-1]:0,r2=i<n-1?v[i+n+1]:0;
    var m11=mix+2*Math.PI*k1*Math.cos(beta*p1[gi]);
    var m22=mix+2*Math.PI*k2*Math.cos(beta*p2[gi]);
    out[i]=(2*u1-l1-r1)*inv+m11*u1+mix*u2;
    out[i+n]=(2*u2-l2-r2)*inv+m22*u2+mix*u1;
  }
}
function hqResidualHalf(p1,p2,g,k1,k2,dx){
  var sp=Math.sqrt(Math.PI), beta=2*sp, Nh=p1.length, n=Nh-2, r=new Float64Array(2*n), mx=0,l2=0;
  for(var i=1;i<Nh-1;i++){
    var j=i-1,d21=(p1[i-1]-2*p1[i]+p1[i+1])/(dx*dx),d22=(p2[i-1]-2*p2[i]+p2[i+1])/(dx*dx),ch=p1[i]+p2[i];
    var r1=-d21+g*g*ch+sp*k1*Math.sin(beta*p1[i]);
    var r2=-d22+g*g*ch+sp*k2*Math.sin(beta*p2[i]);
    r[j]=r1;r[j+n]=r2;mx=Math.max(mx,Math.abs(r1),Math.abs(r2));l2+=r1*r1+r2*r2;
  }
  return{r:r,maxr:mx,l2:Math.sqrt(l2)};
}
function hqCgHalf(b,p1,p2,g,k1,k2,dx,maxIt,tol){
  var N=b.length,x=new Float64Array(N),r=new Float64Array(b),z=new Float64Array(N),p=new Float64Array(N),Ap=new Float64Array(N);
  var Nh=p1.length,n=Nh-2,inv=1/(dx*dx),sp=Math.sqrt(Math.PI),beta=2*sp,mix=g*g;
  function prec(){for(var i=0;i<n;i++){var gi=i+1,d1=2*inv+mix+2*Math.PI*k1*Math.cos(beta*p1[gi]),d2=2*inv+mix+2*Math.PI*k2*Math.cos(beta*p2[gi]);z[i]=r[i]/Math.max(Math.abs(d1),1e-12);z[i+n]=r[i+n]/Math.max(Math.abs(d2),1e-12)}}
  prec();p.set(z);var rz=hqDot(r,z),r0=Math.max(hqNorm(r),1e-300);
  for(var it=0;it<maxIt;it++){
    Ap.fill(0);hqApplyJHalf(p,Ap,p1,p2,g,k1,k2,dx);var den=hqDot(p,Ap);if(!Number.isFinite(den)||Math.abs(den)<1e-300)break;
    var a=rz/den;for(var j=0;j<N;j++){x[j]+=a*p[j];r[j]-=a*Ap[j]}
    if(hqNorm(r)<tol*r0)break;
    prec();var rz2=hqDot(r,z),bcg=rz2/Math.max(rz,1e-300);for(var k=0;k<N;k++)p[k]=z[k]+bcg*p[k];rz=rz2;
  }
  return x;
}
function hqSolveHalfBackground(g,k1,k2,halfL,dx){
  var sp=Math.sqrt(Math.PI), Nh=Math.round(halfL/dx)+1, x=new Array(Nh);for(var i=0;i<Nh;i++)x[i]=-halfL+i*dx;
  dx=x[1]-x[0];var m1=Math.sqrt(2*Math.PI*k1),m2=Math.sqrt(2*Math.PI*k2);
  var p1=x.map(function(xx){return 0.5*sp*(1+Math.tanh(0.5*m1*(xx+2/m1)))});
  var p2=x.map(function(xx){return -0.5*sp*(1+Math.tanh(0.5*m2*(xx+2/m2)))});
  p1[0]=0;p1[Nh-1]=sp/2;p2[0]=0;p2[Nh-1]=-sp/2;
  var dt=0.08*dx*dx;
  for(var it=0;it<250;it++){
    var rr=hqResidualHalf(p1,p2,g,k1,k2,dx).r,n=Nh-2,n1=p1.slice(),n2=p2.slice();
    for(var q=1;q<Nh-1;q++){var j=q-1;n1[q]=p1[q]-dt*rr[j];n2[q]=p2[q]-dt*rr[j+n]}
    p1=n1;p2=n2;p1[0]=0;p1[Nh-1]=sp/2;p2[0]=0;p2[Nh-1]=-sp/2;
  }
  var info=hqResidualHalf(p1,p2,g,k1,k2,dx);
  for(var nw=0;nw<18;nw++){
    if(info.maxr<1e-9)break;
    var b=new Float64Array(info.r.length);for(var bi=0;bi<b.length;bi++)b[bi]=-info.r[bi];
    var del=hqCgHalf(b,p1,p2,g,k1,k2,dx,Math.min(900,6*b.length),1e-9);
    var best1=p1,best2=p2,best=info;
    for(var ls=0;ls<12;ls++){
      var lam=Math.pow(0.5,ls),t1=p1.slice(),t2=p2.slice(),ni=Nh-2;
      for(var qi=1;qi<Nh-1;qi++){var jj=qi-1;t1[qi]=p1[qi]+lam*del[jj];t2[qi]=p2[qi]+lam*del[jj+ni]}
      t1[0]=0;t1[Nh-1]=sp/2;t2[0]=0;t2[Nh-1]=-sp/2;var ti=hqResidualHalf(t1,t2,g,k1,k2,dx);
      if(ti.l2<best.l2){best1=t1;best2=t2;best=ti}if(ti.l2<info.l2*(1-1e-4*lam))break;
    }
    p1=best1;p2=best2;info=best;
  }
  return{x:x,p1:p1,p2:p2,dx:dx,res:info.maxr,resL2:info.l2};
}
function hqMirrorBackground(bg){
  var sp=Math.sqrt(Math.PI),h=bg.x,p1=bg.p1,p2=bg.p2;
  return{x:h.concat(h.slice(0,-1).map(function(v){return -v}).reverse()),p1:p1.concat(p1.slice(0,-1).map(function(v){return sp-v}).reverse()),p2:p2.concat(p2.slice(0,-1).map(function(v){return -sp-v}).reverse()),dx:bg.dx,res:bg.res,resL2:bg.resL2};
}
var oldBuildAtomForPhysicsPatch=buildAtom;
massData=function(p){
  var a=p.g*p.g+2*Math.PI*p.k1,d=p.g*p.g+2*Math.PI*p.k2,b=p.g*p.g,tr=a+d,disc=Math.sqrt((a-d)*(a-d)+4*b*b),lamLight=0.5*(tr-disc),lamHeavy=0.5*(tr+disc),v1=b,v2=lamLight-a,norm=Math.sqrt(v1*v1+v2*v2)||1;v1/=norm;v2/=norm;if(v2<0){v1=-v1;v2=-v2}return{mLight2:Math.max(lamLight,0),mHeavy2:Math.max(lamHeavy,0),mLight:Math.sqrt(Math.max(lamLight,0)),mHeavy:Math.sqrt(Math.max(lamHeavy,0)),vLight:[v1,v2]};
};
potentialGradient=function(phi1,phi2,p){var common=p.g*p.g*(phi1+phi2);return[common+SQRT_PI*p.k1*Math.sin(TWO_SQRT_PI*phi1),common+SQRT_PI*p.k2*Math.sin(TWO_SQRT_PI*phi2)];};
totalEnergy=function(s){var p=s.p,E=0;for(var i=1;i<p.N-1;i++){var d1=(s.phi1[i+1]-s.phi1[i-1])/(2*p.dx),d2=(s.phi2[i+1]-s.phi2[i-1])/(2*p.dx),kin=0.5*(s.pi1[i]*s.pi1[i]+s.pi2[i]*s.pi2[i]),grad=0.5*(d1*d1+d2*d2),pot=0.5*p.g*p.g*(s.phi1[i]+s.phi2[i])*(s.phi1[i]+s.phi2[i])-0.5*p.k1*Math.cos(TWO_SQRT_PI*s.phi1[i])-0.5*p.k2*Math.cos(TWO_SQRT_PI*s.phi2[i]);E+=(kin+grad+pot)*p.dx}return E;};
buildAtom=function(){
  var p=readParams();message('building half-domain Newton atom background...');
  var half=hqSolveHalfBackground(p.g,p.k1,p.k2,0.5*p.L,p.dx);var full=hqMirrorBackground(half);
  p.N=full.x.length;p.dx=full.dx;p.L=full.x[full.x.length-1]-full.x[0];
  var x=new Float64Array(full.x),phi1=new Float64Array(full.p1),phi2=new Float64Array(full.p2),md=massData(p);
  st={p:p,x:x,md:md,phi1Static:phi1.slice(),phi2Static:phi2.slice(),phi1:phi1.slice(),phi2:phi2.slice(),pi1:new Float64Array(p.N),pi2:new Float64Array(p.N),force1:new Float64Array(p.N),force2:new Float64Array(p.N),gamma:makeSponge(p,x),time:0,step:0,running:false,built:true,records:[],initialEnergy:null,lastSpectrum:null,probe:updateProbeIndices(p),backgroundResidual:full.res,backgroundMethod:'half-newton-mirror'};
  if(p.driveMode==='packet')addIncomingPacket(st);st.initialEnergy=totalEnergy(st);updateReadouts(full.res);drawAll();message('atom built by half-domain Newton; residual '+full.res.toExponential(2));
};
ui.buildBtn.removeEventListener('click',oldBuildAtomForPhysicsPatch);
ui.buildBtn.addEventListener('click',buildAtom);
message('ready; click Build atom first, then Solve bound mode');
