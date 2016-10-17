#include "FN_knot_code.h"    //contains some functions and all global variables/*Available options: FROM_B_FILE: Skip initialisation, input from previous run. FROM_K_FILE: Initialise from input file(s) containing points for each link component. TORUS_KNOT: Initialise a p-q torus knot HOPF_LINK: Initialise a Hopf Link BORR_RINGS: Initialise Borromean Rings THREE_TWIST: A 5_2 knot */unsigned int option = TORUS_KNOT;         //unknot default option//const bool addtwist = false;/**If TORUS_KNOT chosen**/int p = 0;int q = 1;/**If FROM_K_FILE chosen**/int ncomp = 1;                          //number of link componentsstring knot_filename = "61_knot";      //assumed input filename format of "XXXXX.txt" for single component or "XXXXXi.txt" for i=1..ncomp for multiple components/**IF FROM_B_FILE chosen**/string B_filename = "phi.vtk";    //filename for B/phi field or uv/**IF HOPF_LINK**/double skew = 1.1;     //Relative link radii, =1 for links of equal size. Choose >=1bool poslink = true;     //True if +1 link no.//Grid pointsunsigned int Nx = 201;   //No. points in x,y and zunsigned int Ny = 201;unsigned int Nz = 141;unsigned int timesteps = 50000;     //total number of timestepsunsigned int skipstep = 2000;       //print out every # timestepsdouble dtime = 0.05;         //size of each time step//Size parameters for knotdouble lambda = 21.3;                //approx wavelengthdouble coresize;   //approx core diameterconst double epsilon = 0.3;                //parameters for F-N eqnsconst double beta = 0.7;const double gam = 0.5;//System size parametersdouble size,h,xmax,ymax,zmax;int NK;   //number of knot points//Unallocated matricesvector<double> X;     //knot position and tangentvector<double> Y;vector<double> Z;vector<double> dlx;vector<double> dly;vector<double> dlz;double length;   //initial knot lengthinline int pt(int i, int j, int k)       //convert i,j,k to single index{    return (i*Ny*Nz+j*Nz+k);}int main (void){    double *x, *y, *z, *phi, *u, *v, *ucv, *d2u, *d2v;    int *missed;    int i,j,k,n,t;        x = new double [Nx];    y = new double [Ny];    z = new double [Nz];    u = new double [Nx*Ny*Nz];    v = new double [Nx*Ny*Nz];    ucv = new double [Nx*Ny*Nz];    d2u = new double [Nx*Ny*Nz];    d2v = new double [Nx*Ny*Nz];    phi = new double [Nx*Ny*Nz];  //scalar potential    missed = new int [Nx*Ny*Nz];        coresize = lambda/(2*M_PI);    size = 10*lambda;    h = size/Nx;    xmax = Nx*h/5;   //max grid point knot can inhabit    ymax = Ny*h/5;   //must be < N*h/2    zmax = Nz*h/5;        for(i=0;i<Nx;i++)           //initialise grid    {        x[i]=h*(i-((int) (Nx-1)/2));    }    for(j=0;j<Ny;j++)    {        y[j]=h*(j-((int) (Ny-1)/2));    }    for(k=0;k<Nz;k++)    {        z[k]=h*(k-((int) (Nz-1)/2));    }        if (option == FROM_B_FILE)    {        cout << "Reading input file...\n";        Bfile_read(phi,missed);    }    else    {        //Initialise knot        length = initialise_knot();        if(length==0)        {            cout << "Error reading input option. Aborting...\n";            return 1;        }                cout << "Total no. knot points: " << NK << endl;        print_knot();                //Calculate B and phi for initial conditions        initial_cond(x,y,z,phi,missed);    }        //Fitzhugh Nagumo equations    cout << "Calculating u and v...\n";    uv_initialise(phi,u,v,ucv,missed);        delete [] missed;    delete [] phi;    delete [] x;    delete [] y;    delete [] z;        double timesec = 60.0*Nx*Ny*Nz/(1000.0*101*101*81*dtime);    unsigned int tt = ((int) timesec);    cout << "Updating u and v: " << tt/60 << "m " << tt%60 << "s per unit T.\n";    for(n=0;n<=timesteps;n++)    {        if(n%skipstep==0)        {            crossgrad_calc(u,v,ucv);            print_uv(u,v,ucv,n*dtime);            cout << "T = " << n*dtime << endl;        }        uv_update(u,v);    }        delete [] u;    delete [] v;    delete [] ucv;    delete [] d2u;    delete [] d2v;        return 0;}/*************************Functions for knot initialisation*****************************/double initialise_knot(){    double L;    switch (option)    {        case FROM_K_FILE: L = init_from_k_file();            break;                    case TORUS_KNOT: L = init_torus_knot();            break;                    case HOPF_LINK:  L = init_hopf_link();            break;                    case BORR_RINGS:  L = init_borr_rings();            break;                    case THREE_TWIST:  L = init_three_twist();            break;           default: L=0;            break;    }        //if (addtwist) L+= add_twist();        return L;}double init_from_k_file(void){    double xt,yt,zt,dx,dy,dz,dl,lseg,Lh;    //temporary variables    vector<double> px,py,pz,dr,ntx,nty,ntz;  //points, distances and tangents    int npts;  //counter    string temp;    double L=0;    int n,m,t,s,NKh;    ifstream knotin;   //knot file(s)    string filename, buff;    stringstream ss;        NK=0;    //count total no. of points        for (m=1; m<=ncomp; m++)    {        npts=0;                ss.clear();        ss.str("");        if (ncomp==1) ss << knot_filename << ".txt";        else ss << knot_filename << m << ".txt";                filename = ss.str();        knotin.open(filename.c_str());                Lh=0;        while(knotin.good())   //read in points for knot        {            if(getline(knotin,buff))            {                ss.clear();                ss.str("");                ss << buff;                ss >> xt >> yt >> zt;            }            else break;                        px.push_back(xt*xmax);            py.push_back(yt*ymax);             //store points            pz.push_back(zt*zmax);                        if(npts>0)            {                dx = px[npts] - px[npts-1];                dy = py[npts] - py[npts-1];     //distance between points                dz = pz[npts] - pz[npts-1];                dr.push_back(sqrt(dx*dx + dy*dy + dz*dz));                ntx.push_back(dx/dr[npts-1]);    //tangent direction to next pt, goes from 0:npts-1                nty.push_back(dy/dr[npts-1]);                ntz.push_back(dz/dr[npts-1]);                Lh += dr[npts-1];                //total length of link component            }            npts++;        }                knotin.close();                NKh = ((int) (2*Lh/h));  //Number of points to define knot with ~h/2 spacing        dl = Lh/NKh;       //Actual spacing ~h/2                //Start at p0        X.push_back(px[0]);        Y.push_back(py[0]);        Z.push_back(pz[0]);                n = 0; //input pt counter                for(t=1;t<NKh;t++)        {            s = NK+t;            X.push_back(X[s-1] + dl*ntx[n]);     //interpolate between input points            Y.push_back(Y[s-1] + dl*nty[n]);            Z.push_back(Z[s-1] + dl*ntz[n]);            lseg = sqrt((X[s] - px[n])*(X[s] - px[n]) + (Y[s] - py[n])*(Y[s] - py[n]) + (Z[s] -               pz[n])*(Z[s] - pz[n]));     //distance from last input point            while(lseg>dr[n])   //if we have passed next input point            {                n=n+1;     //move to next input point                X[s] = px[n] + (lseg-dr[n-1])*ntx[n];    //add extra bit in next direction                Y[s] = py[n] + (lseg-dr[n-1])*nty[n];                Z[s] = pz[n] + (lseg-dr[n-1])*ntz[n];                lseg = sqrt((X[s] - px[n])*(X[s] - px[n]) + (Y[s] - py[n])*(Y[s] - py[n]) + (Z[s] - pz[n])*(Z[s] - pz[n]));    //recalculate segment length            }            //cout << X[s] << ' ' << Y[s] << ' ' << Z[s] << '\n';        }                px.clear();        py.clear();        pz.clear();        ntx.clear();        nty.clear();        ntz.clear();        dr.clear();                //smooth curve        double *Xnew,*Ynew,*Znew;        double d2x,d2y,d2z;                for(n=0;n<10000;n++)        {            for(t=0;t<NKh;t++)            {                d2x = X[NK+incp(t,1,NKh)] - 2*X[NK+t] + X[NK+incp(t,-1,NKh)];                d2y = Y[NK+incp(t,1,NKh)] - 2*Y[NK+t] + Y[NK+incp(t,-1,NKh)];                d2z = Z[NK+incp(t,1,NKh)] - 2*Z[NK+t] + Z[NK+incp(t,-1,NKh)];                X[NK+t] += 0.01*d2x;                Y[NK+t] += 0.01*d2y;                Z[NK+t] += 0.01*d2z;            }        }                for(t=0;t<NKh;t++)        {            dlx.push_back(0.5*(X[NK+incp(t,1,NKh)] - X[NK+incp(t,-1,NKh)]));   //central diff for tangent            dly.push_back(0.5*(Y[NK+incp(t,1,NKh)] - Y[NK+incp(t,-1,NKh)]));            dlz.push_back(0.5*(Z[NK+incp(t,1,NKh)] - Z[NK+incp(t,-1,NKh)]));        }        NK += NKh;        L += Lh;    //keep track of total length and npts to define link    }        return L;}double init_torus_knot(void){    //use xmax as size limit    double L = 0;    double maxlength = 0.5*M_PI*xmax*sqrt(4*q*q+p*p);  //upper bound on knot length    NK = ((int) (5*maxlength/h)); //approximately one knot point every h/5    double dt = 2*M_PI/NK;    //parametric variable    double t;    int n;        for(n=0;n<NK;n++)    {        t=n*dt;        X.push_back(xmax*cos(q*t)*(3+cos(p*t))/4);        Y.push_back(xmax*sin(q*t)*(3+cos(p*t))/4);        Z.push_back(xmax*sin(p*t)/4);        //use analytical expression for dr/dt:        dlx.push_back((-q*sin(q*t)*(3+cos(p*t))-p*sin(p*t)*cos(q*t))*xmax*dt/4);        dly.push_back((q*cos(q*t)*(3+cos(p*t))-p*sin(p*t)*sin(q*t))*xmax*dt/4);        dlz.push_back(p*cos(p*t)*xmax*dt/4);        //don't have analytical expression for total length        L += sqrt((dlx[n])*(dlx[n]) + (dly[n])*(dly[n]) + (dlz[n])*(dlz[n]));    }        return L;}double init_hopf_link(void){    //use zmax as size limit    double r1 = min(0.67*xmax,zmax);    double r2 = min(0.67*xmax,zmax)/skew;    double L1 = 2*M_PI*r1;    double L2 = 2*M_PI*r2;    double x1,y1,z1,x2,y2,z2;    double a1,b1,a2,b2,dt1,dt2,NK1,NK2;    int n,lk;    double t;    double L=L1+L2;    ncomp = 2;        /*First circle*/    x1 = 0.5*r1;    y1 = 0;    z1 = 0;                  //centre coords    NK1 = ((int) 2*L1/h);   //number of points to define circle    dt1 = 2*M_PI/NK1;       //angular spacing of points required        if(poslink) lk=1;    else lk=-1;        /*a1 = 0.01*M_PI;    //arbitrary angles     b1 = 0;          //plane angles with x and y directions     a1 = tan(a1);     b1 = tan(b1);     for(n=0;n<NK1;n++)     {     t=n*dt1;     X.push_back(x1+r1*(-b1*cos(t) - a1*sin(t)/sqrt(1+a1*a1+b1*b1))/sqrt(a1*a1+b1*b1));     Y.push_back(y1+r1*(a1*cos(t) - b1*sin(t)/sqrt(1+a1*a1+b1*b1))/sqrt(a1*a1+b1*b1));     Z.push_back(z1-r1*sqrt(a1*a1+b1*b1)*sin(t)/sqrt(1+a1*a1+b1*b1));     //use analytical expression for dr/dt:     dlx.push_back(dt1*r1*(b1*sin(t) - a1*cos(t)/sqrt(1+a1*a1+b1*b1))/sqrt(a1*a1+b1*b1));     dly.push_back(dt1*r1*(-a1*sin(t) - b1*cos(t)/sqrt(1+a1*a1+b1*b1))/sqrt(a1*a1+b1*b1));     dlz.push_back(-dt1*r1*sqrt(a1*a1+b1*b1)*cos(t)/sqrt(1+a1*a1+b1*b1));     }*/        for(n=0;n<NK1;n++)    {        t=n*dt1;        X.push_back(x1+r1*cos(lk*t));        Y.push_back(y1+r1*sin(lk*t));        Z.push_back(z1);        //use analytical expression for dr/dt:        dlx.push_back(-dt1*lk*r1*sin(lk*t));        dly.push_back(dt1*lk*r1*cos(lk*t));        dlz.push_back(0);    }        /*Second circle*/    x2 = -0.5*r2;    y2 = 0;    z2 = 0;   //centre coords    NK2 = ((int) 2*L2/h);   //number of points to define circle    dt2 = 2*M_PI/NK2;       //angular spacing of points required        /*a2 = 0.01*M_PI;     b2 = 0.51*M_PI;   //plane angles with x and y directions     a2 = tan(a2);     b2 = tan(b2);     for(n=0;n<NK2;n++)     {     t=n*dt2;     X.push_back(x2+r2*(-b2*cos(t) - a2*sin(t)/sqrt(1+a2*a2+b2*b2))/sqrt(a2*a2+b2*b2));     Y.push_back(y2+r2*(a2*cos(t) - b2*sin(t)/sqrt(1+a2*a2+b2*b2))/sqrt(a2*a2+b2*b2));     Z.push_back(z2-r2*sqrt(a2*a2+b2*b2)*sin(t)/sqrt(1+a2*a2+b2*b2));     //use analytical expression for dr/dt:     dlx.push_back(dt2*r2*(b2*sin(t) - a2*cos(t)/sqrt(1+a2*a2+b2*b2))/sqrt(a2*a2+b2*b2));     dly.push_back(dt2*r2*(-a2*sin(t) - b2*cos(t)/sqrt(1+a2*a2+b2*b2))/sqrt(a2*a2+b2*b2));     dlz.push_back(-dt2*r2*sqrt(a2*a2+b2*b2)*cos(t)/sqrt(1+a2*a2+b2*b2));     }*/        for(n=0;n<NK2;n++)    {        t=n*dt2;        X.push_back(x2+r2*cos(t));        Y.push_back(y2);        Z.push_back(z2+r2*sin(t));        //use analytical expression for dr/dt:        dlx.push_back(-dt2*r2*sin(t));        dly.push_back(0);        dlz.push_back(dt2*r2*cos(t));    }        NK = NK1+NK2;        return L;}double init_borr_rings(void)  //not written yet{    return 0;}double init_three_twist(void){    //use xmax as size limit    double L = 0;    double applength = 3*(2*xmax + 3*ymax + 8*zmax);  //approximate knot length    NK = ((int) (5*applength/h)); //approximately one knot point every h/5    double dt = 2*M_PI/NK;    //parametric variable    double t;    int n;        for(n=0;n<NK;n++)    {        t=n*dt;        X.push_back(xmax*cos(2*t+0.2));        Y.push_back(ymax*cos(3*t+0.7));        Z.push_back(zmax*cos(7*t));        //use analytical expression for dr/dt:        dlx.push_back(-2*xmax*sin(2*t+0.2)*dt);        dly.push_back(-3*ymax*sin(3*t+0.7)*dt);        dlz.push_back(-7*zmax*sin(7*t)*dt);        //don't have analytical expression for total length        L += sqrt((dlx[n])*(dlx[n]) + (dly[n])*(dly[n]) + (dlz[n])*(dlz[n]));    }        return L;}/*double add_twist(void)   //only works for a knot at the moment{    double dt=2*M_PI/NK;    double L=0;    double r=0.5*coresize;    double n1[3],n2[3];    int n;    double d2x,d2y,d2z,t;        for(n=0;n<NK;n++)    {        t=n*dt;        d2x=X[incp(n,1,NK)]-2*X[n]+X[incp(n,-1,NK)];        d2y=Y[incp(n,1,NK)]-2*Y[n]+Y[incp(n,-1,NK)];        d2z=Z[incp(n,1,NK)]-2*Z[n]+Z[incp(n,-1,NK)];        n1[0]=d2x/sqrt(d2x*d2x+d2y*d2y+d2z*d2z);        n1[1]=d2y/sqrt(d2x*d2x+d2y*d2y+d2z*d2z);        n1[2]=d2z/sqrt(d2x*d2x+d2y*d2y+d2z*d2z);        n2[0]=(dly[n]*n1[2]-dlz[n]*n1[1])/(dlx[n]*dlx[n] + dly[n]*dly[n] + dlz[n]*dlz[n]);        n2[1]=(dlz[n]*n1[0]-dlx[n]*n1[2])/(dlx[n]*dlx[n] + dly[n]*dly[n] + dlz[n]*dlz[n]);        n2[2]=(dlx[n]*n1[1]-dly[n]*n1[0])/(dlx[n]*dlx[n] + dly[n]*dly[n] + dlz[n]*dlz[n]);        X.push_back(X[n] + r*(cos(10*t)*n1[0] + sin(10*t)*n2[0]));        Y.push_back(Y[n] + r*(cos(10*t)*n1[1] + sin(10*t)*n2[1]));        Z.push_back(Z[n] + r*(cos(10*t)*n1[2] + sin(10*t)*n2[2]));    }        for(n=0;n<NK;n++)    {        dlx.push_back(0.5*(X[NK+incp(n,1,NK)] - X[NK+incp(n,-1,NK)]));   //central diff for tangent        dly.push_back(0.5*(Y[NK+incp(n,1,NK)] - Y[NK+incp(n,-1,NK)]));        dlz.push_back(0.5*(Z[NK+incp(n,1,NK)] - Z[NK+incp(n,-1,NK)]));        L += sqrt((dlx[NK+n])*(dlx[NK+n]) + (dly[NK+n])*(dly[NK+n]) + (dlz[NK+n])*(dlz[NK+n]));    }        NK=2*NK;        return L;}*//*************************Functions for B and Phi calcs*****************************/void initial_cond(double *x, double *y, double *z, double *phi, int *missed){    int *ignore;  //Points to ignore    int *ignore1;    double *Bx;  //Mag field    double *By;    double *Bz;    double *Bmag;    int i,j,k;        ignore = new int [Nx*Ny*Nz];    ignore1 = new int [Nx*Ny*Nz];    Bx = new double [Nx*Ny*Nz];    By = new double [Nx*Ny*Nz];    Bz = new double [Nx*Ny*Nz];    Bmag = new double [Nx*Ny*Nz];        double timesec = 30.0*Nx*Ny*Nz*NK/(101*101*81*1322);    unsigned int tt = ((int) timesec);    cout << "Calculating B field ~ " << tt/60 << "m " << tt%60 << "s...\n";    B_field_calc(x,y,z,Bx, By, Bz, Bmag, ignore, ignore1, missed);    timesec = 130.0*Nx*Ny*Nz*(Nx+Ny+Nz)/(101*101*81*283);    tt = ((int) timesec);    cout << "Calculating scalar potential ~ " << tt/60 << "m " << tt%60 << "s...\n";    phi_calc(Bx, By, Bz, Bmag, ignore, ignore1, missed, phi);    cout << "Printing B and phi...\n";    print_B_phi(x, y, z, Bx, By, Bz, ignore, ignore1, missed, phi);        delete [] ignore;    delete [] ignore1;    delete [] Bx;    delete [] By;    delete [] Bz;    delete [] Bmag;        return;}void B_field_calc(double *x, double *y, double *z, double *Bx, double *By, double *Bz, double *Bmag, int *ignore, int *ignore1, int *missed){    int i,j,k,n,t;    double lx,ly,lz,lmag;        for(i=0;i<Nx;i++)    {        for(j=0;j<Ny;j++)        {            for(k=0;k<Nz;k++)            {                n = pt(i,j,k);    //3D counter                Bx[n] = 0;                By[n] = 0;                Bz[n] = 0;                missed[n] = 1;   //intialise                for(t=0;t<NK;t++)  //integrate over line                {                    lx = x[i]-X[t];    //distance to point on line                    ly = y[j]-Y[t];                    lz = z[k]-Z[t];                    lmag = sqrt(lx*lx + ly*ly + lz*lz);                    if (lmag < 2*coresize) ignore[n]=1;   //do not use these points first time                    if (lmag < 0.5*coresize) ignore1[n]=1; //do not use these at all                    Bx[n] += (ly*dlz[t] - lz*dly[t])/(2*lmag*lmag*lmag);                    By[n] += (lz*dlx[t] - lx*dlz[t])/(2*lmag*lmag*lmag);                    Bz[n] += (lx*dly[t] - ly*dlx[t])/(2*lmag*lmag*lmag);                }                Bmag[n] = sqrt(Bx[n]*Bx[n] + By[n]*By[n] + Bz[n]*Bz[n]);            }        }    }}void phi_calc(double *Bx, double *By, double *Bz, double *Bmag, int *ignore, int *ignore1, int *missed, double *phi){    int i0=(Nx+1)/2;    int j0=(Ny+1)/2;   //base point for path integral    int k0=(Nz+1)/2;    int i[2],j[2],k[2],id,jd,kd,c1,c2,c3,pathlength,t,nt,ntm;    double Bxmid,Bymid,Bzmid;    int *pi,*pj,*pk;    int n = pt(i0,j0,k0);        missed[n]=0;  //matrix to store points where phi is not calculated    phi[n]=0;        pi = new int [Nx+Ny+Nz];    pj = new int [Nx+Ny+Nz];    pk = new int [Nx+Ny+Nz];        for(id=0; id<(Nx+1)/2; id++)    //from zero to half grid points    {        for(jd=0; jd<(Ny+1)/2; jd++)        {            for(kd=0; kd<(Nz+1)/2; kd++)            {                i[0] = id;                i[1] = Nx-1-id;                j[0] = jd;                j[1] = Ny-1-jd;                k[0] = kd;                k[1] = Nz-1-kd;                for(c1=0;c1<2;c1++)    //count inwards from corners                {                    for(c2=0;c2<2;c2++)                    {                        for(c3=0;c3<2;c3++)                        {                            n = pt(i[c1],j[c2],k[c3]);                                                        if(missed[n]==1 && ignore[n]==0)                            {                                pathlength = pathfind(i0,j0,k0,i[c1],j[c2],k[c3],pi,pj,pk,ignore,Bx,By,Bz,Bmag);  //find path to current point                                for (t=1;t<=pathlength;t++)   //travel alog path                                {                                    nt = pt(pi[t],pj[t],pk[t]);     //this point                                    ntm = pt(pi[t-1],pj[t-1],pk[t-1]); //prev pt                                    Bxmid = 0.5*(Bx[nt]+Bx[ntm]);                                    Bymid = 0.5*(By[nt]+By[ntm]);   //midpoint                                    Bzmid = 0.5*(Bz[nt]+Bz[ntm]);                                    phi[nt] = phi[ntm] + h*(Bxmid*(pi[t]-pi[t-1]) + Bymid*(pj[t]-pj[t-1]) + Bzmid*(pk[t]-pk[t-1]));    //integrate along                                    missed[nt]=0;                                    while(phi[nt]>M_PI) phi[nt] -= 2*M_PI;                                    while(phi[nt]<-M_PI) phi[nt] += 2*M_PI;                                }                            }                        }                    }                }            }        }    }        for(id=0; id<Nx; id++)    //fill in ignore points but not ignore1    {        for(jd=0; jd<Ny; jd++)        {            for(kd=0; kd<Nz; kd++)            {                n = pt(id,jd,kd);                if(ignore1[n]==0 && missed[n]==1)                {                    pathlength = pathfind(i0,j0,k0,id,jd,kd,pi,pj,pk,ignore1,Bx,By,Bz,Bmag);                    for (t=1;t<=pathlength;t++)                    {                        nt = pt(pi[t],pj[t],pk[t]);                        ntm = pt(pi[t-1],pj[t-1],pk[t-1]);                        Bxmid = 0.5*(Bx[nt]+Bx[ntm]);                        Bymid = 0.5*(By[nt]+By[ntm]);                        Bzmid = 0.5*(Bz[nt]+Bz[ntm]);                        phi[nt] = phi[ntm] + h*(Bxmid*(pi[t]-pi[t-1]) + Bymid*(pj[t]-pj[t-1]) + Bzmid*(pk[t]-pk[t-1]));                        missed[nt]=0;                        while(phi[nt]>M_PI) phi[nt] -= 2*M_PI;                        while(phi[nt]<-M_PI) phi[nt] += 2*M_PI;                    }                }            }        }    }        delete [] pi;    delete [] pj;    delete [] pk;}int pathfind(int i0, int j0, int k0, int ie, int je, int ke, int *pi, int *pj, int *pk, int *ignore, double *Bx, double *By, double *Bz, double *Bmag){    int io,jo,ko,ip,jp,kp,n,np,nu,go,stop,t=0;    int *track;    double MAX,weight1,weight2;    track = new int [Nx*Ny*Nz];        for(ip=0;ip<Nx;ip++)    {        for(jp=0;jp<Ny;jp++)        {            for(kp=0;kp<Nz;kp++)            {                track[pt(ip,jp,kp)]=0;        //initialise            }        }    }        pi[0] = i0;    //starting point for path    pj[0] = j0;    pk[0] = k0;    int di = ie - i0;    int dj = je - j0;    int dk = ke - k0;   //distance to go to final point        while (t<Nx+Ny+Nz && (abs(di)>0 || abs(dj)>0 || abs(dk)>0))  //until reaches end of path or path too long    {        n = pt(pi[t],pj[t],pk[t]);        nu = pt(pi[t]+sign(di),pj[t]+sign(dj),pk[t]+sign(dk));  //check direct route        if(ignore[nu] + track[nu]==0)  //if next space is available        {            pi[t+1] = pi[t] + sign(di);            pj[t+1] = pj[t] + sign(dj);            pk[t+1] = pk[t] + sign(dk);            t++;   //move to next point            n = pt(pi[t],pj[t],pk[t]);            track[n]=1;        }        else        {            MAX = -10;    //compare point values            go = 0;            for(ip=-1; ip<2; ip++)            {                for(jp=-1; jp<2; jp++)                {                    for(kp=-1; kp<2; kp++)  //check all neighbours                    {                        np = pt(pi[t]+ip,pj[t]+jp,pk[t]+kp);                        if(pi[t]+ip<Nx && pi[t]+ip>0 && pj[t]+jp<Ny && pj[t]+jp>0 && pk[t]+kp<Nz && pk[t]+kp>0) //If it is in the simulation box                        {                            stop = ignore[np] + track[np];  //not allowed to visit ignore points or previously visited points                        }                        else stop = 1;                        if(stop==0)                        {                            go=1;                            //weigting for which point to favour                            //direction of final point weighting                            weight1 = (di*ip + dj*jp + dk*kp)/(sqrt(di*di + dj*dj + dk*dk)*sqrt(ip*ip + jp*jp + kp*kp));                            //direction of B field weighting (helps to choose a direction around a barrier)                            weight2 = (Bx[np]*ip + By[np]*jp + Bz[np]*kp)/(Bmag[np]*sqrt(ip*ip + jp*jp + kp*kp));                            if(weight1 + weight2 > MAX)                            {                                MAX = weight1+weight2;                                io = ip;                                jo = jp;    //store the most favourable point                                ko = kp;                                n = pt(pi[t],pj[t],pk[t]);                                track[n]=1;  //track points visited                            }                        }                    }                }            }            if(go==1)   //found a point to move to            {                pi[t+1] = pi[t]+io;                pj[t+1] = pj[t]+jo;                pk[t+1] = pk[t]+ko;                t++;   //move to next point            }            else            {                if(t==0)                {                    cout << "Could not find path to" << ie << ' ' << je << ' ' << ke << endl;                    return 0;                }                else                {                    t--;  //go back to refind previous point                }            }        }        di = ie - pi[t];        dj = je - pj[t];        dk = ke - pk[t];    }        if (t==Nx+Ny+Nz) t=0; //couldn't find path        delete [] track;        return t;}/*************************Functions for FN dynamics*****************************/void uv_initialise(double *phi, double *u, double *v, double *ucv, int *missed){    int n;        for(n=0; n<Nx*Ny*Nz; n++)    {        u[n] = (1-missed[n])*(2*cos(phi[n]) - 0.4) - missed[n]*0.4;        v[n] = (1-missed[n])*(sin(phi[n]) - 0.4) - missed[n]*0.4;        //if missed set value to -0.4 (vortex centre value)    }}void crossgrad_calc(double *u, double *v, double *ucv){    int i,j,k;    double dxu,dyu,dzu,dxv,dyv,dzv,ucvx,ucvy,ucvz;        for(i=0;i<Nx;i++)    {        for(j=0; j<Ny; j++)        {            for(k=0; k<Nz; k++)   //Central difference            {                dxu = 0.5*(u[pt(incw(i,1,Nx),j,k)]-u[pt(incw(i,-1,Nx),j,k)])/h;                dxv = 0.5*(v[pt(incw(i,1,Nx),j,k)]-v[pt(incw(i,-1,Nx),j,k)])/h;                dyu = 0.5*(u[pt(i,incw(j,1,Ny),k)]-u[pt(i,incw(j,-1,Ny),k)])/h;                dyv = 0.5*(v[pt(i,incw(j,1,Ny),k)]-v[pt(i,incw(j,-1,Ny),k)])/h;                dzu = 0.5*(u[pt(i,j,incw(k,1,Nz))]-u[pt(i,j,incw(k,-1,Nz))])/h;                dzv = 0.5*(v[pt(i,j,incw(k,1,Nz))]-v[pt(i,j,incw(k,-1,Nz))])/h;                ucvx = dyu*dzv - dzu*dyv;                ucvy = dzu*dxv - dxu*dzv;                ucvz = dxu*dyv - dyu*dxv;                ucv[pt(i,j,k)] = sqrt(ucvx*ucvx + ucvy*ucvy + ucvz*ucvz);            }        }    }    }void uv_update(double *u, double *v){    int i,j,k,l,n;    double D2u, **ku, **kv, *uold, *vold;        ku = new double* [4];    kv = new double* [4];    uold = new double [Nx*Ny*Nz];    vold = new double [Nx*Ny*Nz];        for(l=0;l<4;l++)    {        ku[l] = new double [Nx*Ny*Nz];        kv[l] = new double [Nx*Ny*Nz];    }        for(i=0;i<Nx;i++)    {        for(j=0; j<Ny; j++)        {            for(k=0; k<Nz; k++)            {                n = pt(i,j,k);                uold[n] = u[n];  //old value of u                vold[n] = v[n];  //old value of v            }        }    }        for(l=0;l<4;l++)  //u and v update for each fractional time step    {        for(i=0;i<Nx;i++)        {            for(j=0; j<Ny; j++)            {                for(k=0; k<Nz; k++)   //Central difference                {                    n = pt(i,j,k);                    D2u = (u[pt(incw(i,1,Nx),j,k)] + u[pt(incw(i,-1,Nx),j,k)] + u[pt(i,incw(j,1,Ny),k)] + u[pt(i,incw(j,-1,Ny),k)] + u[pt(i,j,incw(k,1,Nz))] + u[pt(i,j,incw(k,-1,Nz))] - 6*u[n])/(h*h);                    ku[l][n] = (u[n] - u[n]*u[n]*u[n]/3 - v[n])/epsilon + D2u;                    kv[l][n] = epsilon*(u[n] + beta - gam*v[n]);                }            }        }                for(i=0;i<Nx;i++)        {            for(j=0; j<Ny; j++)            {                for(k=0; k<Nz; k++)  //update                {                    n = pt(i,j,k);                    if(l==0 || l==1)                    {                        u[n] = uold[n] + 0.5*dtime*ku[l][n];                        v[n] = vold[n] + 0.5*dtime*kv[l][n];                    }                    else                    {                        if(l==2)                        {                            u[n] = uold[n] + dtime*ku[l][n];                            v[n] = vold[n] + dtime*kv[l][n];                        }                        else                        {                            u[n] = uold[n] + dtime*(ku[0][n] + 2*ku[1][n] + 2*ku[2][n] + ku[3][n])/6;                            v[n] = vold[n] + dtime*(kv[0][n] + 2*kv[1][n] + 2*kv[2][n] + kv[3][n])/6;                        }                    }                }            }        }    }        for(l=0;l<4;l++)    {        delete [] ku[l];        delete [] kv[l];    }        delete [] uold;    delete [] vold;    delete [] ku;    delete [] kv;}/*************************File reading and writing*****************************/void print_uv(double *u, double *v, double *ucv, double t){    int i,j,k,n;    stringstream ss;    ss << "uv_plot" << t << ".vtk";    ofstream uvout (ss.str().c_str());        uvout << "# vtk DataFile Version 3.0\nUV fields\nASCII\nDATASET STRUCTURED_POINTS\n";    uvout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';    uvout << "ORIGIN " << -0.5*Nx*h << ' ' << -0.5*Ny*h << ' ' << -0.5*Nz*h << '\n';    uvout << "SPACING " << h << ' ' << h << ' ' << h << '\n';    uvout << "POINT_DATA " << Nx*Ny*Nz << '\n';    uvout << "SCALARS u float\nLOOKUP_TABLE default\n";            for(k=0; k<Nz; k++)    {        for(j=0; j<Ny; j++)        {            for(i=0; i<Nx; i++)            {                n = pt(i,j,k);                uvout << u[n] << '\n';            }        }    }        uvout << "SCALARS ucrossv float\nLOOKUP_TABLE default\n";        for(k=0; k<Nz; k++)    {        for(j=0; j<Ny; j++)        {            for(i=0; i<Nx; i++)            {                n = pt(i,j,k);                uvout << ucv[n] << '\n';            }        }    }        uvout.close();}void print_knot(){    string fn = "knotplot.vtk";    ofstream knotout (fn.c_str());    int t;        knotout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET UNSTRUCTURED_GRID\n";    knotout << "POINTS " << NK << " float\n";    knotout << fixed;    for(t=0;t<NK;t++)    {        knotout << setprecision(5) << X[t] << ' ' << setprecision(5) << Y[t] << ' ' << setprecision(5) << Z[t] << '\n';    }            knotout << "\nPOINT_DATA " << NK << "\nVECTORS Tangent float\n";    for(t=0;t<NK;t++)    {        knotout << setprecision(5) << dlx[t] << ' ' << setprecision(5) << dly[t] << ' ' << setprecision(5) << dlz[t] << '\n';    }        knotout.close();}void print_B_phi(double *x, double *y, double*z, double *Bx, double *By, double *Bz, int *ignore, int *ignore1, int *missed, double *phi){    int i,j,k,n;    string fn = "phi.vtk";        ofstream Bout (fn.c_str());        Bout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET STRUCTURED_POINTS\n";    Bout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';    Bout << "ORIGIN " << -0.5*Nx*h << ' ' << -0.5*Ny*h << ' ' << -0.5*Nz*h << '\n';    Bout << "SPACING " << h << ' ' << h << ' ' << h << '\n';    Bout << "POINT_DATA " << Nx*Ny*Nz << '\n';    Bout << "SCALARS Phi float\nLOOKUP_TABLE default\n";    for(k=0; k<Nz; k++)    {        for(j=0; j<Ny; j++)        {            for(i=0; i<Nx; i++)            {                n = pt(i,j,k);                Bout << phi[n] << '\n';            }        }    }        Bout << "SCALARS Missed float\nLOOKUP_TABLE default\n";    for(k=0; k<Nz; k++)    {        for(j=0; j<Ny; j++)        {            for(i=0; i<Nx; i++)            {                n = pt(i,j,k);                Bout << missed[n] << '\n';            }        }    }        Bout << "VECTORS Bfield float\n";        for(k=0; k<Nz; k++)    {        for(j=0; j<Ny; j++)        {            for(i=0; i<Nx; i++)            {                n = pt(i,j,k);                Bout << Bx[n] << ' ' << By[n] << ' ' << Bz[n] << '\n';            }        }    }        Bout.close();}void Bfile_read(double *phi, int *missed){    int i,j,k,n;    string temp;    ifstream fin (B_filename.c_str());        for(k=0; k<Nz; k++)    {        for(j=0; j<Ny; j++)        {            for(i=0; i<Nx; i++)            {                n = pt(i,j,k);                fin >> temp >> temp >> temp >> temp >> temp >> temp >> phi[n] >> missed[n];            }        }    }        fin.close();}