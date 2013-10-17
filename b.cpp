
double acc(double M,vec r){
     absr=abs(r);
     vec a=G*M/(absr*absr*absr)*r;
     return a;
 }

double RK4(double M,double t,vec z, double dt, dzdt){
    double dtdiv2=dt/2; // dt/2 is used a couple of times
    double thalf=t+dtdiv2; //t+dt/2
    double tnew=t+dt;
    double K1 =acc(M, z);
    double zhalf=z+dtdiv2*K1; //z+dt/2
    double K2=acc(thalf, zhalf);
    double zhalf_=z+dtdiv2*K2; //z+dt/2*K2, _ in zhalf_ instead of z^\tilde
    double K3=acc(thalf,zhalf_);
    double znew=z+K3*dt;
    double K4=acc(tnew,znew);
    znew = z + dt/6.0*(K1+K2+K3+K4);
    return znew;
    }

class solver{


}

int main()
{
    int dim=2; //dimention of problem
    T=5; //
    double G= 6.67428e-11;//gravitational constant
    vec M= zeros(10); //vector for masses
    mat a=zeros(n,);
    r.jupiter()= RK4(M(i),double t,vec z, double dt, dzdt);
    //Sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune, pluto
    M(0)=2.0e30; M(1)=2.4e23; M(2)=4.9e24; M(3)=6.0e24; M(4)=6.6e23;
    M(5)=1.9e27; M(6)=5.5e26; M(7)=8.8e25; M(8)=1.03e26; M(9)=1.31e22;

    return 0;
}
