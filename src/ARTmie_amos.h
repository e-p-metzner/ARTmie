#ifndef ARTMIE_AMOS_HPP
#define ARTMIE_AMOS_HPP

// Real functions

double dgamln(double x);
double dsign(double dest, double src);

// Complex base functions

double zabs(double zr, double zi);
void zmlt(double fac1r, double fac1i, double fac2r, double fac2i, double *OUTR, double *OUTI);
void zdiv(double numr, double numi, double denr, double deni, double *OUTR, double *OUTI);
void zsqrt(double ar, double ai, double *BR, double *BI);
void zexp(double zr, double zi, double *OUTR, double *OUTI);
void zlog(double ar, double ai, double *BR, double *BI, int *IERR);

// Bessel related functions

void zasyi(double zr, double zi, double fnu, int kode, int n, double *YR, double *YI, int *NZ, double rl, double tol, double elim, double alim, int *debug);
void zrati(double zr, double zi, double fnu, int n, double *CYR, double *CYI, double tol, int *debug);
void zshch(double zr, double zi, double *CSHR, double *CSHI, double *CCHR, double *CCHI);
void zs1s2(double zr, double zi, double *S1R, double *S1I, double *S2R, double *S2I, int *NZ, double ascle, double alim, int *IUF, int *debug);
void zuchk(double yr, double yi, int *NZ, double ascle, double tol, int *debug);
void zunhj(double zr, double zi, double fnu, int ipmtr, double tol, double *PHIR, double *PHII, double *ARGR, double *ARGI, \
           double *ZETA1R, double *ZETA1I, double *ZETA2R, double *ZETA2I, double *ASUMR, double *ASUMI, double *BSUMR, double *BSUMI, int *debug);
void zunik(double zr, double zi, double fnu, int ikflg, int ipmtr, double tol, int *INIT, double *PHIR, double *PHII, \
           double *ZETA1R, double *ZETA1I, double *ZETA2R, double *ZETA2I, double *SUMR, double *SUMI, double *CWORKR, double *CWORKI, int *debug);
void zkscl(double zr, double zi, double fnu, int n, double *YR, double *YI, int *NZ, double rzr, double rzi, double ascle, double tol, double elim, int* debug);
void zmlri(double zr, double zi, double fnu, int kode, int n, double *YR, double *YI, int *NZ, double tol, int *debug);
void zseri(double zr, double zi, double fnu, int kode, int n, double *YR, double *YI, int *NZ, double tol, double elim, double alim, int *debug);
void zuoik(double zr, double zi, double fnu, int kode, int ikflg, int n, double *YR, double *YI, int *NUF, double tol, double elim, double alim, int *debug);
void zbknu(double zr, double zi, double fnu, int kode, int n, double *YR, double *YI, int *NZ, double tol, double elim, double alim, int *debug);
void zwrsk(double zr, double zi, double fnu, int kode, int n, double *YR, double *YI, int *NZ, double *CWR, double *CWI, double tol, double elim, double alim, int *debug);
void zacai(double zr, double zi, double fnu, int kode, int mr, int n, double *YR, double *YI, int *NZ, double rl, double tol, double elim, double alim, int *debug);
void zairy(double zr, double zi, int id, int kode, double *AIR, double *AII, int *NZ, int *IERR, int *debug);
void zunk1(double zr, double zi, double fnu, int kode, int mr, int n, double *YR, double *YI, int *NZ, double tol, double elim, double alim, int *debug);
void zunk2(double zr, double zi, double fnu, int kode, int mr, int n, double *YR, double *YI, int *NZ, double tol, double elim, double alim, int *debug);
void zbunk(double zr, double zi, double fnu, int kode, int mr, int n, double *YR, double *YI, int *NZ, double tol, double elim, double alim, int *debug);
void zuni1(double zr, double zi, double fnu, int kode, int n, double *YR, double *YI, int *NZ, int *NLAST, double fnul, double tol, double elim, double alim, int *debug);
void zuni2(double zr, double zi, double fnu, int kode, int n, double *YR, double *YI, int *NZ, int *NLAST, double fnul, double tol, double elim, double alim, int *debug);
void zbuni(double zr, double zi, double fnu, int kode, int n, double *YR, double *YI, int *NZ, int nui, int *NLAST, double fnul, double tol, double elim, double alim, int *debug);
void zbinu(double zr, double zi, double fnu, int kode, int n, double *CYR, double *CYI, int *NZ, double rl, double fnul, double tol, double elim, double alim, int *debug);
void zacon(double zr, double zi, double fnu, int kode, int mr, int n, double *YR, double *YI, int *NZ, double rl, double fnul, double tol, double elim, double alim, int *debug);
void zbesh(double zr, double zi, double fnu, int kode, int m, int n, double *CYR, double *CYI, int *NZ, int *IERR, int *debug);
void zbesi(double zr, double zi, double fnu, int kode, int n, double *CYR, double *CYI, int *NZ, int *IERR, int *debug);
void zbesk(double zr, double zi, double fnu, int kode, int n, double *CYR, double *CYI, int *NZ, int *IERR, int *debug);
void zbesj(double zr, double zi, double fnu, int kode, int n, double *CYR, double *CYI, int *NZ, int *IERR, int *debug);
void zbesy(double zr, double zi, double fnu, int kode, int n, double *CYR, double *CYI, int *NZ, double *CWORKR, double *CWORKI, int *IERR, int *debug);

#endif /* ARTMIE_AMOS_HPP */