#pragma once

#include <random>
#include <stack>
#include <cassert>

#include <sstream>
#include <fstream>

#include <string>
#include <vector>
#include <array>
#include <functional>
// #include <any>


using Idx = std::size_t;
// using Idx = unsigned long int;
// using Idx = int;
using Double = long double;
// using int = short;

using namespace std;

constexpr int TWO = 2;
constexpr int THREE = 3;
constexpr int SIX = 6;

const double abs_tautil = std::sqrt(5.0/2.0); // 1.0 <= abs
const double arg_tautil = std::atan(3.0); // pi/3.0 <= arg <= pi/2.0 // 3.35
double divs[3] = {0.21931939805707432, 0.2500000002279842, 0.25000000000011924};
double div_eps = 0.7627082855646118;


const bool tautil_default = false; // true->below
double tautil1 = 0.3420201433256688;
double tautil2 = 0.9396926207859083 + 0.00001;



int mult = 4; // even
// int mult = 2;
Idx Lx = 3*1*mult; // 12
Idx Ly = 3*1*mult;

// constexpr Idx Lx = 6*4; // 12
// constexpr Idx Ly = 2*Lx;
constexpr int nparallel = 8; // 4; //12


constexpr int nu = 1; // PP, PA, AA, AP
// const Double beta_c = 0.5 * std::log(2.0 + std::sqrt(3.0));
//
// constexpr Double kappa = 2.0/3.0;
// const Double cos6 = std::cos(M_PI/6.0);
// const Double B = cos6 / (1.0 - kappa*kappa * cos6 * cos6);
//
// constexpr Double alat = 1.0/Lx;
// constexpr Double xipsi = std::sqrt( 1.5*std::sqrt(3.0)*alat / (2.0*M_PI) );

#ifndef _OPENMP
// int omp_get_thread_num() { return 0; }
#endif



std::mt19937 gen;
std::uniform_real_distribution<Double> d01D(0.0, 1.0); // (1, 6);
std::uniform_int_distribution<int> d01I(0, 1); // (1, 6);
std::uniform_int_distribution<Idx> d0N(0, Lx*Ly-1); // (1, 6);
void set_gen( const int seed ) {
  std::mt19937 gen0( seed );
  gen.seed( gen0() );
}
Double dist01(){ return d01D(gen); }
Idx dist0N(){ return d0N(gen); }
int distpm1(){ return 2*d01I(gen)-1; }


// ---------------
// GLOBAL FUNCTIONS
// ---------------

double ell0[2]; // = {1.0, 0.0};
double ell1[2]; // = -omega;
double ell2[2]; // -ell2 - ell0

double ell[3]; // {|ell0|, |ell1|, |ell2|}

double ell_star0[2];
double ell_star1[2];
double ell_star2[2];

double e0[2];
double e1[2];
double e2[2];


double kappa[3];
double theta[3];
double cosH[3];
double BetaTilde[3];
double Beta[3];

// double DBetaDMu[3];
// double DBetaDKappa[3];

void set_tautil(){
  if(!tautil_default){
    std::cout << "setting tautil from abs, arg" << std::endl;
    tautil1 = abs_tautil*std::cos(arg_tautil);
    tautil2 = abs_tautil*std::sin(arg_tautil);
  }
}


void set_ell(){
  ell1[0] = 1.0+tautil1;
  ell1[1] = tautil2;

  ell0[0] = 1.0-2.0*tautil1;
  ell0[1] = -2.0*tautil2;

  ell2[0] = -2.0+tautil1;
  ell2[1] = tautil2;

  ell[0] = std::sqrt( ell0[0]*ell0[0] + ell0[1]*ell0[1] );
  ell[1] = std::sqrt( ell1[0]*ell1[0] + ell1[1]*ell1[1] );
  ell[2] = std::sqrt( ell2[0]*ell2[0] + ell2[1]*ell2[1] );
}


void set_ell_star(){
  const double s = 0.5 * (ell[0]+ell[1]+ell[2]);
  const double area = std::sqrt( s*(s-ell[0])*(s-ell[1])*(s-ell[2]) );

  double coeff;
  coeff = 0.25 * (ell[1]*ell[1] + ell[2]*ell[2] - ell[0]*ell[0]) / area;
  ell_star0[0] =  ell0[1] * coeff;
  ell_star0[1] = -ell0[0] * coeff;

  coeff = 0.25 * (ell[2]*ell[2] + ell[0]*ell[0] - ell[1]*ell[1]) / area;
  ell_star1[0] =  ell1[1] * coeff;
  ell_star1[1] = -ell1[0] * coeff;

  coeff = 0.25 * (ell[0]*ell[0] + ell[1]*ell[1] - ell[2]*ell[2]) / area;
  ell_star2[0] =  ell2[1] * coeff;
  ell_star2[1] = -ell2[0] * coeff;
}


void set_kappa(){
  kappa[0] = 2.0*ell[0] / (ell[0] + ell[1] + ell[2]);
  kappa[1] = 2.0*ell[1] / (ell[0] + ell[1] + ell[2]);
  kappa[2] = 2.0*ell[2] / (ell[0] + ell[1] + ell[2]);
}


void set_e(){
  double len;

  len = std::sqrt( ell_star1[0]*ell_star1[0] + ell_star1[1]*ell_star1[1] );
  e1[0] = ell_star1[0]/len;
  e1[1] = ell_star1[1]/len;

  len = std::sqrt( ell_star2[0]*ell_star2[0] + ell_star2[1]*ell_star2[1] );
  e2[0] = ell_star2[0]/len;
  e2[1] = ell_star2[1]/len;

  len = std::sqrt( ell_star0[0]*ell_star0[0] + ell_star0[1]*ell_star0[1] );
  if(len>1.0e-14){
    e0[0] = ell_star0[0]/len;
    e0[1] = ell_star0[1]/len;
  }
  else{
    len = std::sqrt(ell0[0]*ell0[0] + ell0[1]*ell0[1]);
    e0[0] =  ell0[1]/len; // -(e1[0] + e2[0])/std::sqrt(2);
    e0[1] = -ell0[0]/len; // -(e1[1] + e2[1])/std::sqrt(2);
  }
  // len = std::sqrt( ell_star0[0]*ell_star0[0] + ell_star0[1]*ell_star0[1] );
  // e0[0] = ell_star0[0]/len;
  // e0[1] = ell_star0[1]/len;

  // len = std::sqrt( ell_star1[0]*ell_star1[0] + ell_star1[1]*ell_star1[1] );
  // e1[0] = ell_star1[0]/len;
  // e1[1] = ell_star1[1]/len;

  // len = std::sqrt( ell_star2[0]*ell_star2[0] + ell_star2[1]*ell_star2[1] );
  // e2[0] = ell_star2[0]/len;
  // e2[1] = ell_star2[1]/len;
}


void set_theta(){
  theta[0] = std::acos( -e1[0]*e2[0]-e1[1]*e2[1] );
  theta[1] = std::acos( -e2[0]*e0[0]-e2[1]*e0[1] );
  theta[2] = std::acos( -e0[0]*e1[0]-e0[1]*e1[1] );
}


void set_cos(){
  cosH[0] = std::cos(0.5*theta[0]);
  cosH[1] = std::cos(0.5*theta[1]);
  cosH[2] = std::cos(0.5*theta[2]);
}


void set_beta_tilde(){
  BetaTilde[0] = std::atanh( kappa[0]*cosH[1]*cosH[2]/cosH[0] );
  BetaTilde[1] = std::atanh( kappa[1]*cosH[2]*cosH[0]/cosH[1] );
  BetaTilde[2] = std::atanh( kappa[2]*cosH[0]*cosH[1]/cosH[2] );
}

void set_beta(){
  // Beta[0] = -0.5*std::log( std::tanh( BetaTilde[0] ) );
  // Beta[1] = -0.5*std::log( std::tanh( BetaTilde[1] ) );
  // Beta[2] = -0.5*std::log( std::tanh( BetaTilde[2] ) );
  Beta[0] = -0.5*std::log( kappa[0]*cosH[1]*cosH[2]/cosH[0] );
  Beta[1] = -0.5*std::log( kappa[1]*cosH[2]*cosH[0]/cosH[1] );
  Beta[2] = -0.5*std::log( kappa[2]*cosH[0]*cosH[1]/cosH[2] );

  for(int i=0; i<3; i++) if(std::abs(Beta[i])<1.0e-14) Beta[i] = 0.0;
}


// void set_derivs(){
//   for(int i=0; i<3; i++){
//     DBetaDMu[i] = -0.5 * std::tanh(Beta[i]) * std::cosh(Beta[i])*std::cosh(Beta[i]);
//     DBetaDKappa[i] = std::sinh(Beta[i])*std::cosh(Beta[i]) / kappa[i];
//   }
// }


void set_all(){
  set_tautil();
  set_ell();
  set_ell_star();
  set_kappa();
  set_e();
  set_theta();
  set_cos();
  set_beta_tilde();
  set_beta();
  // set_derivs();
}





Idx idx(const Idx x, const Idx y)  { return x + Lx*y; }

void get_xy(Idx& x, Idx& y, const Idx i)  {
  x = (i+Lx)%Lx;
  y = (i-x)/Lx;
}


inline int get_char( const Idx x, const Idx y) { return (x-y+Lx*Ly)%3; } // tri only on char = 1


bool is_site(const Idx x, const Idx y)  {
  const Idx c = get_char(x,y);
  bool res = false;
  if(c==1) res = true; // e.g., (1,0)
  return res;
}

bool is_site(const Idx i)  {
  Idx x,y;
  get_xy(x, y, i);
  return is_site(x,y);
}

bool is_pt_corr(const Idx x, const Idx y)  {
  const Idx c = get_char(x,y);
  bool res = false;
  if(c==0) res = true; // e.g., (1,0)
  return res;
}

bool is_pt_corr(const Idx i)  {
  Idx x,y;
  get_xy(x, y, i);
  return is_pt_corr(x,y);
}


bool is_face(const Idx x, const Idx y)  {
  const Idx c = get_char(x,y);
  bool res = true;
  if(c==1) res = false; // e.g., (1,0)
  return res;
}


bool is_face(const Idx i)  {
  Idx x,y;
  get_xy(x, y, i);
  return is_face(x,y);
}



void cshift_hex(Idx& xp, Idx& yp, const Idx x, const Idx y, const int mu)  {
  int dx = -(mu+2)%3+1;
  int dy = -(mu+1)%3+1;
  if(mu>=3){
    dx *= -1;
    dy *= -1;
  }
  xp = (x+dx+Lx)%Lx;
  yp = (y+dy+Ly)%Ly;
}

void cshift_hex(Idx& ip, const Idx i, const int mu)  {
  Idx x, y, xp, yp;
  get_xy(x, y, i);
  cshift_hex(xp, yp, x, y, mu);
  ip = idx(xp, yp);
}



void cshift_tri(Idx& xp, Idx& yp, const Idx x, const Idx y, const int mu)  {
  int dx, dy;
  if(mu%3==0){
    dx = -1;
    dy = +2;
  }
  else if(mu%3==1){
    dx = -1;
    dy = -1;
  }
  else if(mu%3==2){
    dx = +2;
    dy = -1;
  }
  else assert(false);

  if(mu>=3){
    dx *= -1;
    dy *= -1;
  }
  xp = (x+dx+Lx)%Lx;
  yp = (y+dy+Ly)%Ly;
}


void cshift_tri(Idx& ip, const Idx i, const int mu)  {
  Idx x, y, xp, yp;
  get_xy(x, y, i);
  cshift_tri(xp, yp, x, y, mu);
  ip = idx(xp, yp);
}



// ----------------
// CLASS DEFINITIONS
// ----------------


struct Spin {
  Idx N;
  std::vector<int> s;

  inline int& operator()(const Idx x, const Idx y) { return s[idx(x,y)]; }
  inline int operator()(const Idx x, const Idx y) const { return s[idx(x,y)]; }

  int& operator[](const Idx i) { return s[i]; }
  int operator[](const Idx i) const { return s[i]; }

  Spin() = delete;

  Spin( const int N_ )
    : N(N_)
    , s( N )
  {
    set1();
  }

  void set1() {
    for(Idx i=0; i<Lx*Ly; i++) {
      if(!is_site(i)) s[i]=0;
      else s[i] = 1;
    }
  }

  void random() {
    for(Idx i=0; i<Lx*Ly; i++) {
      if(!is_site(i)) s[i]=0;
      else s[i] = distpm1();
    }
  }


  Double ss_corr( const Idx dx, const Idx dy ) const {
    assert(0<=dx && dx<Lx);
    assert(0<=dy && dy<Ly);

    Double res = 0.0;
    int counter = 0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
	if(!is_site(x,y)) continue;
        const Idx xp = (x+dx)%Lx;
        const Idx yp = (y+dy)%Ly;
	if( !is_site(xp,yp) ) continue;

	res += (*this)(x,y) * (*this)(xp,yp);
	counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> ss_corr() const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
    // #pragma omp parallel for num_threads(nparallel) schedule(static)
#endif
    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
	if(!is_pt_corr(dx, dy)) continue;
        corr[idx(dx,dy)] = ss_corr( dx, dy );
      }}

    return corr;
  }


  Double eps( const Idx x, const Idx y ) const { // defined on face
    assert(0<=x && x<Lx);
    assert(0<=y && y<Ly);
    assert(is_face(x,y));

    Double res = 0.0;

    if(get_char(x,y)==0){
      res += 1 - (*this)( (x-1+Lx)%Lx, (y+1)%Ly ) * (*this)( x, (y-1+Ly)%Ly );
      res += 1 - (*this)( (x+1)%Lx, y ) * (*this)( x, (y-1+Ly)%Ly );
      res += 1 - (*this)( (x+1)%Lx, y ) * (*this)( (x-1+Lx)%Lx, (y+1)%Ly );
    }
    else if(get_char(x,y)==2){
      res += 1 - (*this)( x, (y+1)%Ly ) * (*this)( (x+1)%Lx, (y-1+Ly)%Ly );
      res += 1 - (*this)( x, (y+1)%Ly ) * (*this)( (x-1+Lx)%Lx, y );
      res += 1 - (*this)( (x-1+Lx)%Lx, y ) * (*this)( (x+1)%Lx, (y-1+Ly)%Ly );
    }
    else assert(false);
    res *= 0.25;
    res += -1.0 + div_eps;

    return res;
  }


  Double eps_1pt() const { // both over even/odd
    Double res = 0.0;

// #ifdef _OPENMP
// #pragma omp parallel for collapse(2) num_threads(nparallel)
// #endif
    int counter = 0;
    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_face(x,y) ) continue;
        res += eps( x, y );
        counter++;
      }}

    res /= counter;
    return res;
  }


  Double epseps_corr( const Idx dx, const Idx dy ) const {
    assert(0<=dx && dx<Lx);
    assert(0<=dy && dy<Ly);

    Double res = 0.0;
    int counter = 0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( (x-y+Lx*Ly)%3!=0 ) continue; // from even only
        const Idx xp = (x+dx)%Lx;
        const Idx yp = (y+dy)%Ly;
        if( !is_face(xp,yp) ) continue;

        res += eps(x,y) * eps(xp,yp);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> epseps_corr() const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
    // #pragma omp parallel for num_threads(nparallel) schedule(static)
#endif
    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        if( !is_face(dx,dy) ) continue;
        corr[idx(dx,dy)] = epseps_corr( dx, dy );
      }}

    return corr;
  }


  Double Knew( const int mu_, const Idx x, const Idx y ) const { // both even/odd
    assert(0<=x && x<Lx);
    assert(0<=y && y<Ly);
    assert( is_site(x,y) );

    int mu=mu_;
    // if( !is_link(x,y,mu) ) mu = (mu+3)%6;
    // assert( is_link(x,y,mu) );

    Idx xp, yp;
    cshift_tri( xp, yp, x, y, mu );

    const Double tmp = (*this)(x,y)*(*this)(xp,yp); // @@@
    const Double res = -0.5/kappa[mu%3] * ( tmp - 1.0 ) - divs[mu%3];
    return res;
  }


  Double Knew_1pt( const int mu ) const { // both even/odd
    assert(0<=mu && mu<3);
    // std::vector<Double> tmp(nparallel, 0.0);
    // std::vector<int> counter(nparallel, 0);
    int counter = 0;
    Double res = 0.0;

    // #ifdef _OPENMP
    // #pragma omp parallel for num_threads(nparallel) schedule(static)
    // #endif
    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_site(x,y) ) continue;
        // if( !is_link(x,y,mu) ) continue;
        res += Knew( mu, x, y );
        counter++;
      }}

    // Double res = 0.0;
    // int tot = 0;
    // for(int i=0; i<nparallel; i++) {
    //   res += tmp[i];
    //   tot += counter[i];
    // }
    res /= counter;
    // res /= tot;
    return res;
  }


  Double TM( const int mu, const Idx x, const Idx y ) const { // both even/odd
    assert(0<=x && x<Lx);
    assert(0<=y && y<Ly);
    assert( is_site(x,y) );

    // int mu=mu_;
    // if( !is_link(x,y,mu) ) mu = (mu+3)%6;
    // assert( is_link(x,y,mu) );

    Idx xp, yp;
    Double res = 0.0;

    cshift_hex( xp, yp, x, y, mu );
    res += Knew(mu, x, y);
    if(mu==0) res -= 0.5*( eps( x, (y+1)%Ly ) + eps( (x-1+Lx)%Lx, (y+1)%Ly ) ); // mu deriv
    else if(mu==1) res -= 0.5*( eps( (x-1+Lx)%Lx, y ) + eps( x, (y-1+Ly)%Ly ) ); // mu deriv
    else if(mu==2) res -= 0.5*( eps( (x+1)%Lx, y ) + eps( (x+1)%Lx, (y-1+Ly)%Ly ) ); // mu deriv
    else if(mu==3) res -= 0.5*( eps( x, (y-1+Ly)%Ly ) + eps( (x+1)%Lx, (y-1+Ly)%Ly ) ); // mu deriv
    else if(mu==4) res -= 0.5*( eps( (x+1)%Lx, y ) + eps( x, (y+1)%Ly ) ); // mu deriv
    else if(mu==5) res -= 0.5*( eps( (x-1+Lx)%Lx, y ) + eps( (x-1+Lx)%Lx, (y+1)%Ly ) ); // mu deriv

    // std::cout << "debug. TM = " << res << std::endl;
    // std::cout << "debug. TM(" << x << ", " << y << ")  = " << res << std::endl;
    return res;
  }


  Double Txx( const Idx x, const Idx y ) const {
    assert(0<=x && x<Lx);
    assert(0<=y && y<Ly);

    const Double tanh_x = std::tanh(2.0*Beta[1]);
    const Double tanh_y = std::tanh(2.0*Beta[2]);

    Idx xp, yp;
    Double res = 0.0;
    int mu;

    mu = 1;
    cshift_tri( xp, yp, x, y, mu );
    res += 0.5 * tanh_x  * ( (*this)(x,y)*(*this)(xp,yp) - tanh_x );

    mu = 2;
    cshift_tri( xp, yp, x, y, mu );
    res -= 0.5 * tanh_y  * ( (*this)(x,y)*(*this)(xp,yp) - tanh_y );

    return res;
  }


  Double Txx2( const Idx x, const Idx y ) const {
    assert(0<=x && x<Lx);
    assert(0<=y && y<Ly);

    const Double tanh_x = std::tanh(2.0*Beta[1]);
    const Double tanh_y = std::tanh(2.0*Beta[2]);

    Idx xp, yp;
    Double res = 0.0;
    int mu;

    mu = 1;
    cshift_tri( xp, yp, x, y, mu );
    res += 0.25 * tanh_x  * ( (*this)(x,y)*(*this)(xp,yp) - tanh_x );

    mu = 2;
    cshift_tri( xp, yp, x, y, mu );
    res -= 0.25 * tanh_y  * ( (*this)(x,y)*(*this)(xp,yp) - tanh_y );

    mu = 4;
    cshift_tri( xp, yp, x, y, mu );
    res += 0.25 * tanh_x  * ( (*this)(x,y)*(*this)(xp,yp) - tanh_x );

    mu = 5;
    cshift_tri( xp, yp, x, y, mu );
    res -= 0.25 * tanh_y  * ( (*this)(x,y)*(*this)(xp,yp) - tanh_y );

    return res;
  }

  Double Txy( const Idx x, const Idx y ) const {
    assert(0<=x && x<Lx);
    assert(0<=y && y<Ly);

    const Double tanh_x = std::tanh(2.0*Beta[1]);
    const Double tanh_y = std::tanh(2.0*Beta[2]);

    Idx xp1, yp1;
    Idx xp2, yp2;
    Idx xp1p2, yp1p2;
    Idx xp1m2, yp1m2;
    Idx xm1, ym1;
    Idx xm2, ym2;
    cshift_tri( xp1, yp1, x, y, 1 );
    cshift_tri( xp2, yp2, x, y, 2 );
    cshift_tri( xm1, ym1, x, y, 4 );
    cshift_tri( xm2, ym2, x, y, 5 );
    //
    cshift_tri( xp1p2, yp1p2, xp1, yp1, 2 );
    cshift_tri( xp1m2, yp1m2, xp1, yp1, 5 );

    Double res = 0.0;

    res -= 0.5 * tanh_x*tanh_y  * ( (*this)(x,y)*(*this)(xp1,yp1) - tanh_x ) * ( (*this)(xp1,yp1)*(*this)(xp1p2,yp1p2) - tanh_y );
    res += 0.5 * tanh_x*tanh_y  * ( (*this)(x,y)*(*this)(xp1,yp1) - tanh_x ) * ( (*this)(xp1,yp1)*(*this)(xp1m2,yp1m2) - tanh_y );

    return res;
  }

  Double Txy2( const Idx x, const Idx y ) const {
    assert(0<=x && x<Lx);
    assert(0<=y && y<Ly);

    const Double tanh_x = std::tanh(2.0*Beta[1]);
    const Double tanh_y = std::tanh(2.0*Beta[2]);

    Idx xp1, yp1;
    Idx xp2, yp2;
    Idx xm1, ym1;
    Idx xm2, ym2;
    //
    Idx xp1p2, yp1p2;
    Idx xp1m2, yp1m2;
    //
    Idx xm1p2, ym1p2;
    Idx xm1m2, ym1m2;

    cshift_tri( xp1, yp1, x, y, 1 );
    cshift_tri( xp2, yp2, x, y, 2 );
    cshift_tri( xm1, ym1, x, y, 4 );
    cshift_tri( xm2, ym2, x, y, 5 );
    //
    cshift_tri( xp1p2, yp1p2, xp1, yp1, 2 );
    cshift_tri( xp1m2, yp1m2, xp1, yp1, 5 );
    //
    cshift_tri( xm1p2, ym1p2, xm1, ym1, 2 );
    cshift_tri( xm1m2, ym1m2, xm1, ym1, 5 );

    Double res = 0.0;

    res -= 0.5 * tanh_x*tanh_y  * ( (*this)(x,y)*(*this)(xp1,yp1) - tanh_x ) * ( (*this)(xp1,yp1)*(*this)(xp1p2,yp1p2) - tanh_y );
    res += 0.5 * tanh_x*tanh_y  * ( (*this)(x,y)*(*this)(xp1,yp1) - tanh_x ) * ( (*this)(xp1,yp1)*(*this)(xp1m2,yp1m2) - tanh_y );

    res += 0.5 * tanh_x*tanh_y  * ( (*this)(x,y)*(*this)(xm1,ym1) - tanh_x ) * ( (*this)(xm1,ym1)*(*this)(xm1p2,ym1p2) - tanh_y );
    res -= 0.5 * tanh_x*tanh_y  * ( (*this)(x,y)*(*this)(xm1,ym1) - tanh_x ) * ( (*this)(xm1,ym1)*(*this)(xm1m2,ym1m2) - tanh_y );

    res -= 0.5 * tanh_x*tanh_y  * ( (*this)(x,y)*(*this)(xp2,yp2) - tanh_y ) * ( (*this)(xp2,yp2)*(*this)(xp1p2,yp1p2) - tanh_x );
    res += 0.5 * tanh_x*tanh_y  * ( (*this)(x,y)*(*this)(xp2,yp2) - tanh_y ) * ( (*this)(xp2,yp2)*(*this)(xm1p2,ym1p2) - tanh_x );

    res += 0.5 * tanh_x*tanh_y  * ( (*this)(x,y)*(*this)(xm2,ym2) - tanh_y ) * ( (*this)(xm2,ym2)*(*this)(xp1m2,yp1m2) - tanh_x );
    res -= 0.5 * tanh_x*tanh_y  * ( (*this)(x,y)*(*this)(xm2,ym2) - tanh_y ) * ( (*this)(xm2,ym2)*(*this)(xm1m2,ym1m2) - tanh_x );

    res /= 4.0;

    return res;
  }


  Double TM_1pt( const int mu ) const { // both even/odd
    assert(0<=mu && mu<3);
    // std::vector<Double> tmp(nparallel, 0.0);
    // std::vector<int> counter(nparallel, 0);
    int counter = 0;
    Double res = 0.0;

    // #ifdef _OPENMP
    // #pragma omp parallel for num_threads(nparallel) schedule(static)
    // #endif
    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_site(x,y) ) continue;
        res += TM( mu, x, y );
        counter++;
        // tmp[omp_get_thread_num()] += K( x, y, mu );
        // counter[omp_get_thread_num()]++;
      }}

    // Double res = 0.0;
    // int tot = 0;
    // for(int i=0; i<nparallel; i++) {
    //   res += tmp[i];
    //   tot += counter[i];
    // }
    res /= counter;
    // res /= tot;
    return res;
  }


  Double Txx_1pt() const { // both even/odd
    int counter = 0;
    Double res = 0.0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_site(x,y) ) continue;
        res += Txx( x, y );
        counter++;
      }}

    res /= counter;
    return res;
  }


  Double Txx2_1pt() const { // both even/odd
    int counter = 0;
    Double res = 0.0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_site(x,y) ) continue;
        res += Txx2( x, y );
        counter++;
      }}

    res /= counter;
    return res;
  }

  Double Txy_1pt() const { // both even/odd
    int counter = 0;
    Double res = 0.0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_site(x,y) ) continue;
        res += Txy( x, y );
        counter++;
      }}

    res /= counter;
    return res;
  }


  Double Txy2_1pt() const { // both even/odd
    int counter = 0;
    Double res = 0.0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_site(x,y) ) continue;
        res += Txy2( x, y );
        counter++;
      }}

    res /= counter;
    return res;
  }


  Double TMTN_corr( const int mu, const int nu,
                    const Idx dx, const Idx dy ) const { // only from even
    assert(0<=dx && dx<Lx);
    assert(0<=dy && dy<Ly);
    assert(0<=mu && mu<3);
    assert(0<=nu && nu<3);

    Double res = 0.0;
    int counter = 0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_site(x,y) ) continue;
        // if( (x-y+Lx*Ly)%3!=0 ) continue;
        // assert( is_link(x,y,mu) );
        //

        const Idx xp = (x+dx)%Lx;
        const Idx yp = (y+dy)%Ly;

        // if( !is_link(xp,yp,nu) ) continue;
        if( !is_site(xp,yp) ) continue;
        res += TM(mu,x,y) * TM(nu,xp,yp);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> TMTN_corr(const int mu, const int nu) const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
    // #pragma omp parallel for num_threads(nparallel) schedule(static)
#endif
    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        // if( !is_site(dx,dy) ) continue;
        if(!is_pt_corr(dx, dy)) continue;
        // if( (dx-dy+Lx*Ly)%3!=0 ) continue;
        corr[idx(dx,dy)] = TMTN_corr( mu, nu, dx, dy );
      }}

    return corr;
  }


  Double TxxTxx_corr( const Idx dx, const Idx dy ) const { // only from even
    assert(0<=dx && dx<Lx);
    assert(0<=dy && dy<Ly);

    Double res = 0.0;
    int counter = 0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_site(x,y) ) continue;

        const Idx xp = (x+dx)%Lx;
        const Idx yp = (y+dy)%Ly;

        if( !is_site(xp,yp) ) continue;
        res += Txx(x,y) * Txx(xp,yp);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> TxxTxx_corr() const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
#endif
    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        if(!is_pt_corr(dx, dy)) continue;
        corr[idx(dx,dy)] = TxxTxx_corr( dx, dy );
      }}

    return corr;
  }


  Double Txx2Txx2_corr( const Idx dx, const Idx dy ) const { // only from even
    assert(0<=dx && dx<Lx);
    assert(0<=dy && dy<Ly);

    Double res = 0.0;
    int counter = 0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_site(x,y) ) continue;

        const Idx xp = (x+dx)%Lx;
        const Idx yp = (y+dy)%Ly;

        if( !is_site(xp,yp) ) continue;
        res += Txx2(x,y) * Txx2(xp,yp);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> Txx2Txx2_corr() const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
#endif
    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        if(!is_pt_corr(dx, dy)) continue;
        corr[idx(dx,dy)] = Txx2Txx2_corr( dx, dy );
      }}

    return corr;
  }



  Double TxyTxy_corr( const Idx dx, const Idx dy ) const { // only from even
    assert(0<=dx && dx<Lx);
    assert(0<=dy && dy<Ly);

    Double res = 0.0;
    int counter = 0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_site(x,y) ) continue;

        const Idx xp = (x+dx)%Lx;
        const Idx yp = (y+dy)%Ly;

        if( !is_site(xp,yp) ) continue;
        res += Txy(x,y) * Txy(xp,yp);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> TxyTxy_corr() const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
#endif
    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        if(!is_pt_corr(dx, dy)) continue;
        corr[idx(dx,dy)] = TxyTxy_corr( dx, dy );
      }}

    return corr;
  }


  Double Txy2Txy2_corr( const Idx dx, const Idx dy ) const { // only from even
    assert(0<=dx && dx<Lx);
    assert(0<=dy && dy<Ly);

    Double res = 0.0;
    int counter = 0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_site(x,y) ) continue;

        const Idx xp = (x+dx)%Lx;
        const Idx yp = (y+dy)%Ly;

        if( !is_site(xp,yp) ) continue;
        res += Txy2(x,y) * Txy2(xp,yp);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> Txy2Txy2_corr() const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
#endif
    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        if(!is_pt_corr(dx, dy)) continue;
        corr[idx(dx,dy)] = Txy2Txy2_corr( dx, dy );
      }}

    return corr;
  }


  Double TxxTxy_corr( const Idx dx, const Idx dy ) const { // only from even
    assert(0<=dx && dx<Lx);
    assert(0<=dy && dy<Ly);

    Double res = 0.0;
    int counter = 0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_site(x,y) ) continue;

        const Idx xp = (x+dx)%Lx;
        const Idx yp = (y+dy)%Ly;

        if( !is_site(xp,yp) ) continue;
        res += Txx(x,y) * Txy(xp,yp);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> TxxTxy_corr() const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
#endif
    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        if(!is_pt_corr(dx, dy)) continue;
        corr[idx(dx,dy)] = TxxTxy_corr( dx, dy );
      }}

    return corr;
  }


  Double Txx2Txy2_corr( const Idx dx, const Idx dy ) const { // only from even
    assert(0<=dx && dx<Lx);
    assert(0<=dy && dy<Ly);

    Double res = 0.0;
    int counter = 0;

    for(Idx x=0; x<Lx; x++){
      for(Idx y=0; y<Ly; y++){
        if( !is_site(x,y) ) continue;

        const Idx xp = (x+dx)%Lx;
        const Idx yp = (y+dy)%Ly;

        if( !is_site(xp,yp) ) continue;
        res += Txx2(x,y) * Txy2(xp,yp);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> Txx2Txy2_corr() const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
#endif
    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        if(!is_pt_corr(dx, dy)) continue;
        corr[idx(dx,dy)] = Txx2Txy2_corr( dx, dy );
      }}

    return corr;
  }




  Double TM_ss( const int mu,
                const Idx x0, const Idx y0,
                const Idx x1, const Idx y1,
                const Idx x2, const Idx y2 ) const {
    assert(0<=mu && mu<3);
    assert(0<=x0 && x0<Lx);
    assert(0<=y0 && y0<Ly);
    assert(0<=x1 && x1<Lx);
    assert(0<=y1 && y1<Ly);
    assert(0<=x2 && x2<Lx);
    assert(0<=y2 && y2<Ly);

    Double res = 0.0;
    int counter = 0;

    // for(Idx x=0; x<Lx; x++){
    //   for(Idx y=0; y<Ly; y++){
    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        const Idx x0p = (dx+x0)%Lx;
        const Idx y0p = (dy+y0)%Ly;
        const Idx x1p = (dx+x1)%Lx;
        const Idx y1p = (dy+y1)%Ly;
        const Idx x2p = (dx+x2)%Lx;
        const Idx y2p = (dy+y2)%Ly;

        // int mu=mu_;
        // const int c=get_char(x0p,y0p);
        // if(c==2) mu+=3;
        // else if(c==1) continue;
        // assert( is_link(x0p,y0p,mu) );

        // if( !is_link(x0p,y0p,mu) ) continue;
        if( !is_site(x0p,y0p) || !is_site(x1p,y1p) || !is_site(x2p,y2p) ) continue;
        res += TM(mu,x0p,y0p) * (*this)(x1p,y1p) * (*this)(x2p,y2p);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> TM_ss_corr( const int mu,
                                  const Idx x1, const Idx y1, const Idx x2, const Idx y2 ) const {
    assert(0<=mu && mu<3);
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
    // #pragma omp parallel for num_threads(nparallel) schedule(static)
#endif
    for(Idx x0=0; x0<Lx; x0++){
      for(Idx y0=0; y0<Ly; y0++){
        if( !is_site(x0,y0) ) continue;
        corr[idx(x0,y0)] = TM_ss( mu, x0, y0, x1, y1, x2, y2 );
      }}
    return corr;
  }





  Double Txx_ss( const Idx x0, const Idx y0,
                 const Idx x1, const Idx y1,
                 const Idx x2, const Idx y2 ) const {
    assert(0<=x0 && x0<Lx);
    assert(0<=y0 && y0<Ly);
    assert(0<=x1 && x1<Lx);
    assert(0<=y1 && y1<Ly);
    assert(0<=x2 && x2<Lx);
    assert(0<=y2 && y2<Ly);

    Double res = 0.0;
    int counter = 0;

    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        const Idx x0p = (dx+x0)%Lx;
        const Idx y0p = (dy+y0)%Ly;
        const Idx x1p = (dx+x1)%Lx;
        const Idx y1p = (dy+y1)%Ly;
        const Idx x2p = (dx+x2)%Lx;
        const Idx y2p = (dy+y2)%Ly;

        if( !is_site(x0p,y0p) || !is_site(x1p,y1p) || !is_site(x2p,y2p) ) continue;
        res += Txx(x0p,y0p) * (*this)(x1p,y1p) * (*this)(x2p,y2p);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> Txx_ss_corr( const Idx x1, const Idx y1, const Idx x2, const Idx y2 ) const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
    // #pragma omp parallel for num_threads(nparallel) schedule(static)
#endif
    for(Idx x0=0; x0<Lx; x0++){
      for(Idx y0=0; y0<Ly; y0++){
        if( !is_site(x0,y0) ) continue;
        corr[idx(x0,y0)] = Txx_ss( x0, y0, x1, y1, x2, y2 );
      }}
    return corr;
  }


  Double Txx2_ss( const Idx x0, const Idx y0,
                  const Idx x1, const Idx y1,
                  const Idx x2, const Idx y2 ) const {
    assert(0<=x0 && x0<Lx);
    assert(0<=y0 && y0<Ly);
    assert(0<=x1 && x1<Lx);
    assert(0<=y1 && y1<Ly);
    assert(0<=x2 && x2<Lx);
    assert(0<=y2 && y2<Ly);

    Double res = 0.0;
    int counter = 0;

    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        const Idx x0p = (dx+x0)%Lx;
        const Idx y0p = (dy+y0)%Ly;
        const Idx x1p = (dx+x1)%Lx;
        const Idx y1p = (dy+y1)%Ly;
        const Idx x2p = (dx+x2)%Lx;
        const Idx y2p = (dy+y2)%Ly;

        if( !is_site(x0p,y0p) || !is_site(x1p,y1p) || !is_site(x2p,y2p) ) continue;
        res += Txx2(x0p,y0p) * (*this)(x1p,y1p) * (*this)(x2p,y2p);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> Txx2_ss_corr( const Idx x1, const Idx y1, const Idx x2, const Idx y2 ) const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
    // #pragma omp parallel for num_threads(nparallel) schedule(static)
#endif
    for(Idx x0=0; x0<Lx; x0++){
      for(Idx y0=0; y0<Ly; y0++){
        if( !is_site(x0,y0) ) continue;
        corr[idx(x0,y0)] = Txx2_ss( x0, y0, x1, y1, x2, y2 );
      }}
    return corr;
  }



  Double Txy_ss( const Idx x0, const Idx y0,
                 const Idx x1, const Idx y1,
                 const Idx x2, const Idx y2 ) const {
    assert(0<=x0 && x0<Lx);
    assert(0<=y0 && y0<Ly);
    assert(0<=x1 && x1<Lx);
    assert(0<=y1 && y1<Ly);
    assert(0<=x2 && x2<Lx);
    assert(0<=y2 && y2<Ly);

    Double res = 0.0;
    int counter = 0;

    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        const Idx x0p = (dx+x0)%Lx;
        const Idx y0p = (dy+y0)%Ly;
        const Idx x1p = (dx+x1)%Lx;
        const Idx y1p = (dy+y1)%Ly;
        const Idx x2p = (dx+x2)%Lx;
        const Idx y2p = (dy+y2)%Ly;

        if( !is_site(x0p,y0p) || !is_site(x1p,y1p) || !is_site(x2p,y2p) ) continue;
        res += Txy(x0p,y0p) * (*this)(x1p,y1p) * (*this)(x2p,y2p);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> Txy_ss_corr( const Idx x1, const Idx y1, const Idx x2, const Idx y2 ) const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
    // #pragma omp parallel for num_threads(nparallel) schedule(static)
#endif
    for(Idx x0=0; x0<Lx; x0++){
      for(Idx y0=0; y0<Ly; y0++){
        if( !is_site(x0,y0) ) continue;
        corr[idx(x0,y0)] = Txy_ss( x0, y0, x1, y1, x2, y2 );
      }}
    return corr;
  }


  Double Txy2_ss( const Idx x0, const Idx y0,
                  const Idx x1, const Idx y1,
                  const Idx x2, const Idx y2 ) const {
    assert(0<=x0 && x0<Lx);
    assert(0<=y0 && y0<Ly);
    assert(0<=x1 && x1<Lx);
    assert(0<=y1 && y1<Ly);
    assert(0<=x2 && x2<Lx);
    assert(0<=y2 && y2<Ly);

    Double res = 0.0;
    int counter = 0;

    for(Idx dx=0; dx<Lx; dx++){
      for(Idx dy=0; dy<Ly; dy++){
        const Idx x0p = (dx+x0)%Lx;
        const Idx y0p = (dy+y0)%Ly;
        const Idx x1p = (dx+x1)%Lx;
        const Idx y1p = (dy+y1)%Ly;
        const Idx x2p = (dx+x2)%Lx;
        const Idx y2p = (dy+y2)%Ly;

        if( !is_site(x0p,y0p) || !is_site(x1p,y1p) || !is_site(x2p,y2p) ) continue;
        res += Txy2(x0p,y0p) * (*this)(x1p,y1p) * (*this)(x2p,y2p);
        counter++;
      }}

    res /= counter;
    return res;
  }


  std::vector<Double> Txy2_ss_corr( const Idx x1, const Idx y1, const Idx x2, const Idx y2 ) const {
    std::vector<Double> corr(N, 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
    // #pragma omp parallel for num_threads(nparallel) schedule(static)
#endif
    for(Idx x0=0; x0<Lx; x0++){
      for(Idx y0=0; y0<Ly; y0++){
        if( !is_site(x0,y0) ) continue;
        corr[idx(x0,y0)] = Txy2_ss( x0, y0, x1, y1, x2, y2 );
      }}
    return corr;
  }







  std::string print() const {
    std::stringstream ss;
    // for(Idx y=Ly-1; y>=0; y--){
    for(Idx y=0; y<Ly; y++){// y=Ly-1; y>=0; y--){
      for(Idx x=0; x<Lx; x++) {
        ss << std::setw(5) << (*this)(x, y);
      }
      ss << std::endl;
    }
    return ss.str();
  }



  void ckpoint( const std::string& str ) const {
    std::ofstream of( str, std::ios::out | std::ios::binary | std::ios::trunc);
    if(!of) assert(false);

    int tmp = 0.0;
    for(Idx i=0; i<Lx*Ly; i++){
      tmp = (*this)[i];
      of.write( (char*) &tmp, sizeof(int) );
    }
    of.close();
  }

  void read( const std::string& str ) {
    std::ifstream ifs( str, std::ios::in | std::ios::binary );
    if(!ifs) assert(false);

    int tmp;
    for(Idx i=0; i<Lx*Ly; ++i){
      ifs.read((char*) &tmp, sizeof(int) );
      (*this)[i] = tmp;
    }
    ifs.close();
  }



};


void heatbath( Spin& s ){
  // omp not allowed
  for(Idx i=0; i<Lx*Ly; i++){
    if( !is_site(i) ) continue;

    double henv = 0;
    for(int mu=0; mu<SIX; mu++){
      Idx j;
      cshift_tri( j, i, mu );
      henv += Beta[mu%THREE]*s[j];
    }

    const Double p = std::exp(2.0*henv);
    const Double r = dist01();
    if( r<p/(1.0+p) ) s[i] = 1;
    else s[i] = -1;
  }
}


void wolff( Spin& s ){
  std::vector<int> is_cluster(Lx*Ly, 0);
  std::stack<Idx> stack_idx;

  Idx init = dist0N();
  while( !is_site(init) ) init = dist0N();

  is_cluster[init] = 1;
  stack_idx.push(init);

  while( stack_idx.size() != 0 ){

    const Idx p = stack_idx.top();
    stack_idx.pop();
    s[p] = -s[p]; // flip when visited

    for(int mu = 0; mu < SIX; mu++){
      // if( !is_link(p,mu) ) continue;
      Idx q;
      cshift_tri(q, p, mu);
      if( s[q] == s[p] || is_cluster[q]==1 ) continue; // s[x]*sR[y]<0 or y in c

      const Double r = dist01();
      if( r < std::exp(-2.0 * Beta[mu%THREE]) ) continue; // reject

      is_cluster[q] = 1;
      stack_idx.push(q);
    }
  }
}



struct Scalar {
  Double v;

  Scalar()
    : v(0.0)
  {}

  Scalar( const Double v_ )
    : v(v_)
  {}

  Scalar( const Scalar& other )
    : v(other.v)
  {}

  void clear(){ v = 0.0; }

  Scalar& operator+=(const Scalar& rhs){
    v += rhs.v;
    return *this;
  }

  Scalar& operator+=(const Double& rhs){
    v += rhs;
    return *this;
  }

  Scalar& operator/=(const Double& rhs){
    v /= rhs;
    return *this;
  }

  std::string print(const std::function<bool(const Idx x, const Idx y)> is_measure) const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(15);
    ss << v;
    return ss.str();
  }

  void print(std::FILE* stream,
	     const std::function<bool(const Idx x, const Idx y)> is_measure) const {
    fprintf( stream, "%0.15Le\t", v );
  }



};



struct Corr {
  std::vector<Double> v;
  // std::function<bool(const Idx x, const Idx y)> is_measure;

  Corr()
    : v(Lx*Ly, 0.0)
    // , is_measure(is_measure_)
  {}

  Corr( const std::vector<Double> v_ )
    : v(v_)
    // , is_measure(is_measure_)
  {}

  Corr( const Corr& other )
    : v(other.v)
    // , is_measure(other.is_measure)
  {}

  Double& operator()(const Idx x, const Idx y) { return v[idx(x,y)]; }
  Double operator()(const Idx x, const Idx y) const { return v[idx(x,y)]; }

  void clear(){ for(Idx i=0; i<Lx*Ly; i++) v[i] = 0.0; }

  Corr& operator+=(const Corr& rhs)
  {
    assert( rhs.v.size()==Lx*Ly );
    for(Idx i=0; i<Lx*Ly; i++) v[i] += rhs.v[i];
    return *this;
  }

  Corr& operator*=(const double& a)
  {
    for(Idx i=0; i<Lx*Ly; i++) v[i] *= a;
    return *this;
  }

  Corr& operator+=(const double& a)
  {
    for(Idx i=0; i<Lx*Ly; i++) v[i] += a;
    return *this;
  }


  Corr& operator+=(const std::vector<Double>& rhs)
  {
    assert( rhs.size()==Lx*Ly );
    for(Idx i=0; i<Lx*Ly; i++) v[i] += rhs[i];
    return *this;
  }

  Corr& operator/=(const Double& rhs)
  {
    for(Idx i=0; i<Lx*Ly; i++) v[i] /= rhs;
    return *this;
  }


  // void print() const {
  //   for(int y=0; y<Ly; y++){
  //     for(int x=0; x<Lx; x++) {
  //       printf( "%0.15e\t", (*this)(x, y) );
  //     }
  //     printf("\n");
  //   }
  // }

  void print(std::FILE* stream,
	     const std::function<bool(const Idx x, const Idx y)> is_measure ) const {
    for(Idx y=0; y<Ly; y++){
      for(Idx x=0; x<Lx; x++) {
	// if( !is_site(x,y) ) continue;
	// if( get_char(x,y)!=0 ) continue;
	if(!is_measure(x, y)) continue;
        fprintf( stream, "%0.15Le\t", (*this)(x, y) );
      }
      fprintf( stream, "\n");
    }
  }

  std::string print( const std::function<bool(const Idx x, const Idx y)> is_measure ) const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(15);
    for(Idx y=0; y<Ly; y++){
      for(Idx x=0; x<Lx; x++) {
	// if( !is_site(x,y) ) continue;
	// if( get_char(x,y)!=0 ) continue;
	if(!is_measure(x, y)) continue;
        ss << (*this)(x, y) << " ";
      }
      ss << std::endl;
    }
    return ss.str();
  }

};


template<typename T> // T needs to have: .clear, +=, /= defined
struct Obs {
  std::string description;
  int N;
  std::function<T(const Spin&)> f;

  T sum;
  int counter;

  std::function<bool(const Idx x, const Idx y)> is_measure;

  // Obs() = delete;

  Obs
  (
   const std::string& description_,
   const int N_,
   const std::function<T(const Spin&)>& f_,
   const std::function<bool(const Idx x, const Idx y)> is_measure_=[](const Idx x, const Idx y){ return true; }
   )
    : description(description_)
    , N(N_)
    , f(f_)
    , sum()
    , counter(0)
    , is_measure(is_measure_)
  {}

  void clear(){
    sum.clear();
    counter = 0;
  }

  void meas( const Spin& s ) {
    sum += f(s);
    counter++;
  }

  // T mean() const {
  //   T tmp(sum);
  //   tmp /= counter;
  //   return mean;
  // }

  void write_and_clear( const std::string& dir, const int label ){
    const std::string filename = dir + description + "_" + std::to_string(label) + ".dat";

    // std::ofstream of( filename, std::ios::out | std::ios::trunc );
    // if(!of) assert(false);
    // of << std::scientific << std::setprecision(15);
    // of << sum.print();
    // of.close();

    FILE *stream = fopen(filename.c_str(), "w");
    if (stream == NULL) assert(false);
    std::ofstream of( filename );
    sum /= counter;
    sum.print( stream, is_measure);
    fclose( stream );

    clear();
  }

};




