double zernike (float rho,float theta,int n,int m,
		float *dz_drho,float *dz_dtheta);
double zernike_distortions (float rho,float theta,float *coeffs,int nz,
			    float *dz_drho_tot,float *dz_dtheta_tot);

#ifdef ZERNIKE_HOME
int z_edge_amplitudes=0;
#else
extern int z_edge_amplitudes;
#endif
