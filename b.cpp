#include <iostream>
#include <armadillo>
#include <vector>
using namespace std;
using namespace arma;

// http://ssd.jpl.nasa.gov/?constants:
double AU =149.6e6;  // [km] i en astronomisk lengde-enhet
double day = 86400; // [s] 'Julian day'
double Msun= 1988500e24; // [kg]
double G =6.67e-11/( (AU*1.0e3)*(AU*1.0e3)*(AU*1.0e3) )*day*day*Msun; // [m^3 s^-2 kg^-1] =>
// Data for posisjoner og hastigheter fra NASA er gitt i [km].
// Da er det greit aa ha AU i [km] videre i programmet.

class planeter{
    public:
    double r0x;
    double r0y;
    double v0x;
    double v0y;
    double m;
    void set_values(double, double, double,double,double);
};

void planeter::set_values(double R0X, double R0Y,double V0X, double V0Y, double M){
    // verdiene under faar SI-enheter hvis posisjon, hastighet og masse gis i
    // [km], [km/s] og [kg].
    r0x=R0X/AU;//[m]
    r0y=R0Y/AU;
    v0x=V0X/AU*day;//[m/s]
    v0y=V0Y/AU*day;
    m=M/Msun;
}

class systemet{
    public:
    int N;
    int k;
    int task;
    double dt;
    double dtdiv2;
    double dtdiv6;
    double total_mass;
    double L;
    double E_k;
    double E_pot;
    vec m;
    mat r;
    mat v;
    void konstanter(int,int,int,double,double,double,vec);
    void posvel(mat,mat);
    mat doubleRK4();
    mat acc(mat);
    mat massesenter();
    vec energi();
    void bevarte_stoerrelser();
};

void systemet::konstanter(int dimensjoner, int ant_planeter, int tall,double step, double stephalf, double stepsixth, vec M){
    N=dimensjoner; // x, y (, z)
    k=ant_planeter;
    dt=step; // tidssteg
    dtdiv2= stephalf; // dt/2
    dtdiv6=stepsixth; // dt/6
    m=M; // vektor som innehaalder massene på legemene (inkludert sola)
    total_mass=0; // summen av massene til legemene
    for(int i=0;i<k;i++){
        total_mass += m(i);
    }
    task=tall;
}

void systemet::posvel(mat position, mat velocity){
    // Brukes bare til aa initialisere r og v;
    r= position; // Posisjonsmatrise: Hver rad er posisjonsvektoren til et legeme.
    v= velocity; // Hastighetsmatrise: Hver rad er hastighetsvektoren til et legeme.
}

// aksellerasjonsfunksjonen acc(posisjon): regner ut en aksellerasjonsmatrise (aksellerasjonen til legemene)
mat systemet::acc(mat posisjon){
    mat a=zeros<mat>(k,N); //aksellerasjonsmatrise: Hver rad vil inneholde aksellerasjonsvektorer til et av legemene.
                           // faareloepig satt til en k x N nullmatrise
    for(int i=0;i<k;i++){

        for(int j=i+1;j<k;j++){

        // to loekker for aa loepe gjennom alle planetkombinasjoner.
        // j>i slik at vi:
        //       (1): ikke tar med samme planetpar to ganger
        //       (2): (i!=j) slik at legemene vi velger blir forskjellige
            // Tar vektordifferansen mellom de to legemene vi valgte ut
            mat r_vec = posisjon.row(i)-posisjon.row(j); // avstandsvektoren peker på legeme i fra legeme j
            double r_abs = sqrt( r_vec(0)*r_vec(0)+r_vec(1)*r_vec(1) ); // |r_ij|=sqrt(x^2+y^2)
            double forceconstant=G/(r_abs*r_abs*r_abs); // k= G/|r|^3
            double forceconstanti=-m(j)*forceconstant; // Disse 'konstantene' regnes ut paa forhaand
            double forceconstantj=m(i)*forceconstant;  // for spare operasjoner lengre ned.
            if( (i!=1) && (task==0 ) ){
                forceconstanti=0*forceconstant;
            }


            // Merk at fortegnet settes i disse to linjene!
            // a_ij skal peke på legeme j. r_vec peker på planet i. Da må a(i,:) vaere antiparalell
            // med r_vec og a(j,:) være paralell med r_vec.
            // Dette er ordnet i forceconstanti og forceconstantj.

            // Newtons 2. lov gir at kraften mellom de to legemene er like store men motsatt rettet.
            // Aksellerasjonene er også motsatt rettet, men man må gange med massen til legemet
            // man regner ut aksellerasjonen til(se 5,6 linjer over).

            a.row(i)+= forceconstanti*r_vec;  // Summerer aksellerasjonen siden det bare er en faktor m_i som skiller
            a.row(j)+= forceconstantj*r_vec;  // mellom kraftsummen og aksellerasjonssummen. sum_j( F_ij ) = m_i* sum_i(a_i)

        }

    }
    //cout <<"109:  f="<< a.row(1)*m(1)<<endl;
    return a;
}
mat systemet::doubleRK4(){
    // integratormetoden Runge-Kutta 4.
    // funksjonen 'doubleRK4' integrerer hastighet og posisjon paralellt(dobbeltintegral).
    // derav 'double' i navnet på funksjonen

    mat V1= acc(r);
    mat R1= v;
    mat vhalf=v+dtdiv2*V1;
    mat rhalf=r+dtdiv2*R1;

    mat V2= acc(rhalf);
    mat R2= vhalf;
    mat vhalf_=v+dtdiv2*V2;
    mat rhalf_=r+dtdiv2*R2;

    mat V3= acc(rhalf_);
    mat R3= vhalf_;
    mat vnew=v+dt*V3;
    mat rnew=r+dt*R3;

    mat V4= acc(rnew);
    mat R4= vnew;
    v=v+dtdiv6*(V1+2*V2+2*V3+V4);
    r=r+dtdiv6*(R1+2*R2+2*R3+R4);//endring
    if(task!=1){//hvis solen er fiksert i origo
        v.row(0)=zeros<mat>(1,N);
        r.row(0)=zeros<mat>(1,N);
    }
    return r;
}

// 'massesenter(...)' regner ut massesenteret i solsystemet.
mat systemet::massesenter(){
    mat CM= zeros<mat>(1,N);

    for(int i=0;i<k;i++){
        CM+= m(i)*r.row(i);
        }

    CM/=total_mass;
    return CM;
}

// 'bevarte_stoerrelser()' regner ut energi og angulaermoment for maale
// hvor stabil metoden vaar er:
void systemet::bevarte_stoerrelser(){
    // ------------ initialisering ------------
    E_k=0; // total kinetisk energi
    E_pot=0; //total potensiell energi (0 ved uendelig avstand)
    L = 0; //totalt angulaer-moment
    // ----------------------------------------

    for(int i=0;i<k;i++){ //for hver planet:

        //--------------- regner ut E_k og  L  --------------
        double v2 = 0; // |v|^2, kvadratet av v

        //finn |v|^2:
        for(int dim=0;dim<N;dim++){
            // summerer kvadratene av hver hastighetskomponent:
            v2 += v(i,dim)*v(i,dim); //|v|^2 = (delta x)^2 + (delta y)^2 (evt. + (delta z)^2 )
        }

        L += m(i)*(r(i,0)*v(i,1)-r(i,1)*v(i,0)); //Antar at massesenteret er konstant: CM=(0,0)
                                                 // Forenkler legemene til punktlegemer i utregningen av L
        // legger til kinetisk energi for planet i
        // til summen av kinetisk energi til alle legemene:
        E_k += m(i)*v2; // E_k(i) =1/2 * m_i * v^2
        //-------------------------------------------

        //---------------- regner ut E_pot -----------
        for(int j=i+1;j<k;j++){
            // finner par av planeter (bare en gang for hvert par!):
            // avstandsvektor r_vec:
            mat r_vec=r.row(i)-r.row(j);
            // skalar avstand:
            double r_abs = sqrt( r_vec(0)*r_vec(0)+r_vec(1)*r_vec(1) );// (delta x)^2 + (delta y)^2 (evt. + (delta z)^2 )
            // |r_i-r_j|, avstand mellom planet i og j
            E_pot -= G*m(i)*m(j)/r_abs; // Potensiell energi er null ved uendelig avstand
        }
        //------------------------------------------------
    }
    E_k/=2;
}

int main(){
    cout.precision(16);

    cout << "skal sola stå i ro? ja: 0, nei: 1"<<endl <<"for oppgave d: 2"<<endl;
    int task;
    cin >> task;

    //Hvis man skal legge til ytteligere objekter (for eksempel 'rocket'), må man legge til:
    // rocket i listen under og legge til:
    // rocket.set_values(x [km], y [km], v0x[km/s], v0y[km/s], mass);
    // Saa maa man legge til enda en iftest under lengre ned i programmet.

    // objekter i klassen initialverdier: Alle planeter legges til selv om vi ikke bruker alle.

    planeter sun, mercury,venus, earth, mars, jupiter, saturn, uranus, neptune, pluto,rocket;
    sun.set_values(0.0, 0.0, 0.0, 0.0, 1988500e24);// Data jeg hentet fra NASA var relativt til solen.
    mercury.set_values(4.186328367979202E+07,-4.639488984492734E+07,2.649100893125037E+01, 3.496360017130548E+01,0.3302e24);
    venus.set_values(9.107058018777086E+07, -5.942962959851662E+07, 1.892621227368753E+01, 2.918642007512648E+01,48.685e23);

    if(task==0){
        earth.set_values( 1.49597870700E+08, 0, 0, 29.772, 5.97219e24);//oppgave c
    }
    if(task==1){// standard
        earth.set_values(1.401314535653816E+08, 5.135065473668917E+07, -1.074300988397089E+01, 2.785075899706422E+01, 5.97219e24);
    }
    if(task==2){
        earth.set_values( 1.49597870700E+08, 0, 0, 42.1, 5.97219e24);//oppgave d
    }
    mars.set_values(-1.224534489688664E+08, 2.113985089132638E+08, -2.004968649521087E+01, -1.008362569300063E+01 ,6.4185e23);
    jupiter.set_values(-1.100325291552387E+08, 7.644745914375829E+08, -1.310155701837505E+01, -1.242878877322166E+00, 1898.13e24);//original m=1898.13e24
    saturn.set_values(-1.072870309112834E+09, -1.011008964947510E+09, 6.090460191409890E+00, -7.056500798998084E+00,5.68319e26);
    uranus.set_values(2.948078131100535E+09, 5.430764934436175E+08, -1.293465578556454E+00, 6.377544079287039E+00,86.8103e24);
    neptune.set_values(4.032757532991157E+09, -1.962430340115208E+09, 2.332592485566963E+00, 4.917320763804554E+00,102.41e24);
    pluto.set_values(8.986889310434605E+08, -4.775385695669049E+09, 5.439588543356710E+00, -1.008189815772441E-01, 1.309e22);
    rocket.set_values(2.401314535653816E+08, 6.135065473668917E+07, 3.074300988397089E+01, 2.85075899706422E+01, 1.0e6);
    int n;
    cout<<"skriv inn antall punkter i arrayen, int n="<<endl;
    cin >> n;

    double T;
    cout<<"skriv inn maks antall dager for simuleringen"<<endl;
    cin >> T;

    cout<<"skriv inn antall legemer(inkludert sola) i systemet, int k="<<endl;
    int k;
    cin>> k;
    vec tallrekke =zeros(k); // Vektor med tall tillhørede legemene vi velger aa ta med. Solen (nr. 0) er alltid inkludert.
    cout<< "Velg hvilke planeter som skal vaere med. Skriv inn tallet tilhoerende planetene i stigende rekkefoelge."<<endl<<endl;
    cout<< "Mercury: 1" << endl << "venus: 2"<<endl << "earth: 3"<< endl<< "mars: 4"<<endl<<"jupiter: 5"<<endl;
    cout<< "saturn: 6"<<endl <<"uranus: 7"<<endl<<"neptune: 8"<< endl <<"pluto: 9"<<endl<<"rocket: 10"<<endl;

    for(int i=1; i<k;i++){
        cout <<"tall paa neste planet du vil ha med:";
        int a;
        cin>>a;
        tallrekke(i)=a;
    }

    int N=2;//dimensjoner
    double dt=T/(n-1);
    vec t=linspace(0,T,n);
    double dtdiv2=dt/2;//brukes i Runge Kutta
    double dtdiv6=dt/6;// -||-

    vec M=zeros(k); //vektor som inneholder massen til legemene
    mat V= mat(k,N);// hastighetsmatrise med legemenes hastighetsvektorer som rader
    mat R= mat(k,N);// posisjonmatrise med legemenes posisjonsvektorer som rader

    const char *a[k];

    a[0] = "Sun";
    M(0)=sun.m;
    V(0,0)=sun.v0x; V(0,1)=sun.v0y;
    R(0,0)=sun.r0x; R(0,1)=sun.r0y;

    for(int i=1;i<k;i++){
        int tall= tallrekke(i);
        if (tall==1){
            M(i)=mercury.m;
            V(i,0)=mercury.v0x; V(i,1)=mercury.v0y;
            R(i,0)=mercury.r0x; R(i,1)=mercury.r0y;
            a[i]="Mercury";
        }
        if (tall==2){
            M(i)=venus.m;
            V(i,0)=venus.v0x; V(i,1)=venus.v0y;
            R(i,0)=venus.r0x; R(i,1)=venus.r0y;
            a[i] = "Venus";
        }
        if (tall==3){
            M(i)=earth.m;
            V(i,0)=earth.v0x; V(i,1)=earth.v0y;
            R(i,0)=earth.r0x; R(i,1)=earth.r0y;
            a[i]="Earth";
        }
        if (tall==4){
            M(i)=mars.m;
            V(i,0)=mars.v0x; V(i,1)=mars.v0y;
            R(i,0)=mars.r0x; R(i,1)=mars.r0y;
            a[i]="Mars";
        }
        if (tall==5){
            M(i)=jupiter.m;
            V(i,0)=jupiter.v0x; V(i,1)=jupiter.v0y;
            R(i,0)=jupiter.r0x; R(i,1)=jupiter.r0y;
            a[i]="Jupiter";
        }
        if (tall==6){
            M(i)=saturn.m;
            V(i,0)=saturn.v0x; V(i,1)=saturn.v0y;
            R(i,0)=saturn.r0x; R(i,1)=saturn.r0y;
            a[i]="Saturn";
        }
        if (tall==7){
            M(i)=uranus.m;
            V(i,0)=uranus.v0x; V(i,1)=uranus.v0y;
            R(i,0)=uranus.r0x; R(i,1)=uranus.r0y;
            a[i]="Uranus";
        }
        if (tall==8){
            M(i)=neptune.m;
            V(i,0)=neptune.v0x; V(i,1)=neptune.v0y;
            R(i,0)=neptune.r0x; R(i,1)=neptune.r0y;
            a[i]="Neptune";
        }
        if (tall==9){
            M(i)=pluto.m;
            V(i,0)=pluto.v0x; V(i,1)=pluto.v0y;
            R(i,0)=pluto.r0x; R(i,1)=pluto.r0y;
            a[i]="Pluto";
        }
//#if 0
        if (tall==10){
            M(i)=rocket.m;
            V(i,0)=rocket.v0x; V(i,1)=rocket.v0y;
            R(i,0)=rocket.r0x; R(i,1)=rocket.r0y;
            a[i]="rocket";
        }
//#endif
    }

    systemet solar; // Lager et objekt kalt solar
    // konstanter tilhoerende systemet.
    solar.konstanter(N, k, task, dt, dtdiv2, dtdiv6, M);
    solar.posvel(R, V);
    if (task==1){
        //--------------------------------------------------------
        // Justerer hastighetene med hastigheten til massesenteret.
        // Da blir det mye enklere å regne ut angulaermoment.

        mat CM=solar.massesenter();
        mat V_CM= mat(1,N);
        for(int i=0;i<k;i++){
        V_CM += solar.m(i)/(solar.total_mass) * V.row(i);
        }
        for(int i=0;i<k;i++){
            R.row(i)-=CM;
            V.row(i)-=V_CM; // trekker fra hastighetsvektoren til
                           // massesenteret fra hver hastighetsvektor
        }
        solar.posvel(R, V); // fra naa av heter det solar.r og solar.v
        //------------------------------------------------------------
    }
    int delta;
    cout << "hvor mange dt skal det vaere mellom hver gang programmet skriver ut?";
    cin>> delta;
    int counter=delta;// skal printe ut for t=0. Se if-test under

    string filnavn="data.txt";
    string filnavn2="data2.txt";
    ofstream myfile, yourfile;
    myfile.open (filnavn.c_str());
    yourfile.open (filnavn2.c_str());

    // --- print konstanter til fil ---------
    myfile << k <<"  "<< n <<"  "<< N << "  " << T<<endl;
    for(int i=0;i<k;i++){
        myfile << a[i]<<" ";
    }
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    myfile << endl;
    mat printvector= zeros<mat>(1,1+N*k); // en for t, N*k for posisjon
    solar.bevarte_stoerrelser();
    //yourfile << solar.E_k<< "  "<< solar.E_pot<< "  " << solar.L<<endl;
    for(int i =0;i<n;i++){
        if(counter==delta){ //printer ut etter delta dager
            printvector(0)=t(i);
            for(int j=0;j<k;j++){
                printvector(N*j+1)=solar.r(j,0);
                printvector(N*j+2)=solar.r(j,1);
            }

            myfile << printvector; // faar alt på en linje
            solar.bevarte_stoerrelser();
            yourfile << solar.E_k<<"  "<< solar.E_pot <<"  "<< solar.L<<endl;
            counter=0;
        }
        counter+=1;

        solar.doubleRK4();// integrer ett steg
    }

    // gjoer meg sikker på at jeg printer ut data for siste datapunkt:
    printvector(0)=t(n-1);
    for(int j=0;j<k;j++){
        printvector(N*j+1)=solar.r(j,0);
        printvector(N*j+2)=solar.r(j,1);
    }
    myfile << printvector; // faar alt på en linje
    solar.bevarte_stoerrelser();
    yourfile << solar.E_k<<"  "<< solar.E_pot <<"  "<< solar.L;
    // ----------------------------------------------------------------
    myfile.close();
    yourfile.close();
    return 0;
}
