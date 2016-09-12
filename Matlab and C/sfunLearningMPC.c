#define S_FUNCTION_NAME  sfunLearningMPC
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>

/* NPSOL function prototypes */
#include "f2c.h"
#include "npsol.h"

/* MPC and vehicle parameters, utilities */
#include "data_vehicle.h"
#include "data_mpc.h"
#include "mpc_utils.h"

/***** DECLARATIONS *****/
/*** Global variables ***/
static double ZOlGlobal[nz*(Hp+1)];
static double ZRefGlobal[nz*(Hp+1)], URefGlobal[nu*Hu];
static double zCurrGlobal[nz], uPrevGlobal[nu];
static double costZ, derivCost, costU, costdU, costLane;
static double costSlack, costZSplit[nz], costZTerm, costZTermCons[nz-1];
static double Q[nz], R[nu], P[nu], T[nSlack];
static double bLGlobal[NCTOTL], bUGlobal[NCTOTL];
static double coeffTermCost[nPolyOrderTermCost+1];
static double coeffTermCons[(nz-1)][(nPolyOrderTermCons+1)];
static double coeffTermCost_1[nPolyOrderTermCost+1];
static double coeffTermCons_1[(nz-1)][(nPolyOrderTermCons+1)];
static double epsTermConsTol[(nz-1)];
static double eY_LB, eY_UB;
static double coeffRoadAngle[nPolyOrderRoadAngle+1];
static double coeffRoadCurv[nPolyOrderRoadCurv+1];
static double coeffVy[4];
static double coeffYawRate[3];
static double coeffVx[3];
static int trial_number;
static int ReachedSurface;
static double s_start;
static int useQuadSlackCost;
static time_T t_curr;

/*** Inputs ***/
enum { INP_ENABLE, INP_STATE, INP_U_CURR, INP_TUNING, INP_U_LIN, INP_U_REF,
              INP_Z_REF, INP_COEFF_TERM_COST, INP_COEFF_TERM_CONS,
              INP_COEFF_ROAD_ANGLE, INP_COEFF_ROAD_CURV, INP_TRIAL_NUMBER, INP_S_START, INP_COEFF_Vy, 
              INP_COEFF_YawRate, INP_COEFF_Vx, n_Inputs };

/*** Outputs ***/
enum { OUT_DELTA, OUT_AX, OUT_STAT_OPT, OUT_SLACK, OUT_PARAMS, OUT_U_LIN_S,
       OUT_Z_OL, OUT_STATE ,n_Outputs };

/*** PWork vectors ***/
enum { n_PWorks };

/***** MPC SUBROUTINES *****/
/*** Continuous time dynamic bicycle model ***/
/* Point mass models
 * States: [xDot, yDot, PsiDot, Psi, Y, X]'
 * Inputs: beta
 * Parameters: V
 */
static void dynamicBicycleCoeff(double *ZNext, double *z, double *u)
{
    int k;
    double alpha_f, alpha_r, Fc_f, Fc_r, dsdt;

    /* states and inputs */
    double xDot = z[0];
    double yDot = z[1];
    double psiDot = z[2];
    double ePsi = z[3];
    double eY = z[4];
    double s = z[5];
    double af = z[6];
    double delta_f = u[0];
    double ax = u[1];

    /* road angle */
    //double ePsi = psi - polyval(coeffRoadAngle, nPolyOrderRoadAngle, s);

    /* road curvature */
    double road_curv = polyval(coeffRoadCurv, nPolyOrderRoadCurv, s);

    /* Slip angles */
    alpha_f = -delta_f + atan((yDot + l_f*psiDot)/max_value(1.0,xDot));
    alpha_r = atan((yDot - l_r*psiDot)/max_value(1.0,xDot));

    /* Cornering forces (linear) */
    Fc_f = -C_alpha_f*alpha_f;
    Fc_r = -C_alpha_r*alpha_r;

    /* State derivatives */
    dsdt = (xDot*cos(ePsi) - yDot*sin(ePsi))/(1-eY*road_curv);
    //ZNext[0] = xDot   + (yDot*psiDot + af) * dt;
    ZNext[0] = coeffVx[2]*xDot   + coeffVx[0]*yDot + coeffVx[1]*psiDot + af;
    ZNext[1] = yDot + coeffVy[0] * yDot / xDot + 
             + coeffVy[1] * psiDot * xDot + coeffVy[2] * psiDot / xDot + coeffVy[3] * delta_f; 
    //ZNext[1] = yDot   + ( -xDot*psiDot + (2/m)*(Fc_f*cos(delta_f) + Fc_r) ) * dt;
    ZNext[2] = psiDot + coeffYawRate[0] * psiDot / xDot +
             + coeffYawRate[1] * yDot / xDot  + coeffYawRate[2] * delta_f;
    //ZNext[2] = psiDot + ( (2/Iz)*(l_f*Fc_f - l_r*Fc_r) ) * dt;
    ZNext[3] = ePsi   + ( psiDot - dsdt*road_curv ) * dt;
    ZNext[4] =   eY   + ( xDot*sin(ePsi) + yDot*cos(ePsi) ) * dt;
    ZNext[5] =    s   + ( dsdt ) * dt;
//    ZNext[6] =   af + (- af/0.8   + ax/0.8)*dt;
    ZNext[6] =   af + (- af/0.1   + ax/0.1)*dt;
    //Delay
}

static void dynamicBicycle(double *Zdot, double *z, double *u)
{
    int k;
    double alpha_f, alpha_r, Fc_f, Fc_r, dsdt;

    /* states and inputs */
    double xDot = z[0];
    double yDot = z[1];
    double psiDot = z[2];
    double ePsi = z[3];
    double eY = z[4];
    double s = z[5];
    double af = z[6];
    double delta_f = u[0];
    double ax = u[1];

    /* road angle */
    //double ePsi = psi - polyval(coeffRoadAngle, nPolyOrderRoadAngle, s);

    /* road curvature */
    double road_curv = polyval(coeffRoadCurv, nPolyOrderRoadCurv, s);

    /* Slip angles */
    alpha_f = -delta_f + atan((yDot + l_f*psiDot)/max_value(1.0,xDot));
    alpha_r = atan((yDot - l_r*psiDot)/max_value(1.0,xDot));

    /* Cornering forces (linear) */
    Fc_f = -C_alpha_f*alpha_f;
    Fc_r = -C_alpha_r*alpha_r;

    /* State derivatives */
    dsdt = (xDot*cos(ePsi) - yDot*sin(ePsi))/(1-eY*road_curv);
    Zdot[0] = yDot*psiDot + af;
    Zdot[1] = -xDot*psiDot + (2/m)*(Fc_f*cos(delta_f) + Fc_r);
    Zdot[2] = (2/Iz)*(l_f*Fc_f - l_r*Fc_r);
    Zdot[3] = psiDot - dsdt*road_curv;
    Zdot[4] = xDot*sin(ePsi) + yDot*cos(ePsi);
    Zdot[5] = dsdt;
    Zdot[6] = ax/0.1 - af/0.1;
}

/*** Discrete time kinematic bicycle model ***/
/* Point mass models
 * States: [s, eY, Psi, V]'
 * Inputs: [delta, ax]
 */

static void kinematicBicycle(double *zNext, double *z, double *u)
{
    int k;
    double dsdt;

    /* states and inputs */
    double s = z[0];
    double eY = z[1];
    double ePsi = z[2];
    double xDot = z[3];
    double beta_f = atan(l_r/(l_f+l_r)*tan(u[0])); // CHECK!
    double ax = u[1];

    /* road angle */
    //double ePsi = psi - polyval(coeffRoadAngle, nPolyOrderRoadAngle, s);

    /* road curvature */
    double road_curv = polyval(coeffRoadCurv, nPolyOrderRoadCurv, s);

    /* State derivatives */
    dsdt = xDot*cos(ePsi+beta_f)/(1-eY*road_curv);
    zNext[0] = z[0] + dt*dsdt;
    zNext[1] = z[1] + dt*xDot*sin(ePsi+beta_f);
    zNext[2] = z[2] + dt*( xDot/l_r*sin(beta_f) - dsdt*road_curv);
    zNext[3] = z[3] + dt*ax;
}

/*** Simulate vehicle trajectory using point-mass longitudinal model ***/
static void simBicycleModel(double *Zcurr, double *U)
{
    int j, k;
    double t, ZNext[nz], z[nz], Zdot[nz];

    /* Initialization */
    for (k = 0; k < nz; k++){
        ZOlGlobal[k] = Zcurr[k];
    }

    /* Propagation of dynamics */
    for (j = 0; j < Hp; j++){
        /* Get current state */
        for (k = 0; k < nz; k++){
            z[k] = ZOlGlobal[j*nz+k];
        }

        /* Low speed check */
        // if (z[0] < 1.0){
        //     z[0] = 1.0;
        // }

        /* First-order Euler integration */
        
        if ((trial_number == 1)){
        t = 0;
        while (t < dt){
            dynamicBicycle(Zdot, z, U + my_floor(((double) j)/((double) Hi))*nu);
            for (k = 0; k < nz; k++){
                z[k] += dtInt*Zdot[k];
            }
            t = t + dtInt;
        }
            for (k = 0; k < nz; k++){
                ZOlGlobal[(j+1)*nz+k] = z[k];
            }
        }else{
            dynamicBicycleCoeff(ZNext, z, U + my_floor(((double) j)/((double) Hi))*nu);
            for (k = 0; k < nz; k++){
                ZOlGlobal[(j+1)*nz+k] = ZNext[k];
            }
        }
        
        /*
        t = 0;
        while (t < dt){
            dynamicBicycle(Zdot, z, U + my_floor(((double) j)/((double) Hi))*nu);
            for (k = 0; k < nz; k++){
                z[k] += dtInt*Zdot[k];
            }
            t = t + dtInt;
        }
        for (k = 0; k < nz; k++){
                ZOlGlobal[(j+1)*nz+k] = z[k];
        }
        */
    }
}

/***** NPSOL Subroutines *****/
/*** Function implementing the objective function for NPSOL ***/
int *funobj(int *MODE, int *N, double *X, double *OBJ, double *OBJGRAD, int *NSTATE)
{
    int i,j;
    double U[nHu], currCostZ, ParInt;
    double termStateErr;
    double DsDt, derivCost;
    
    /* Initialization */
    // X contains the optimizaton variables, which are the control inputs (a, d_f) and the ParInt-value for parameter interpolation in the Q-function   
    for (i = 0; i < nHu; i++){
        U[i] = X[i];
    }
    ParInt = X[nHu];
    
    /* Simulation */
    simBicycleModel(zCurrGlobal,U);

    /* Objective function value */
    *OBJ = 0;
    /* State costs and lane constraints */
    derivCost = costZ = costLane = 0;
    init_zeros(costZSplit, nz);
    for (i = 0; i < Hp; i++){       // for all steps i=1 to N
        for (j = 0; j < nz; j++){   // and for all states j=1:4
            /* Stage cost */
            //if ( ((trial_number > 1) ||( ZOlGlobal[(i+1)*nz+5]+ s_start - ZRefGlobal[(i+1)*nz+5] > 0 ) )){
            // The following are derivative costs (costs on fast changing variables)
            // it guarantees that the solution is smooth (Ugos Master thesis?)
            if ( ((trial_number > 1) )){        // if at least lap number 2
                if ((j==0)&&( i<(Hp-1) )){      // if state 0 (=s) and pred.<N-1
                    //was 100
                    derivCost = 100 * pow( (ZOlGlobal[(i+1)*nz+j] - ZOlGlobal[(i+1+1)*nz+j] ),2);
                }
                if (( j==1 )&&( i<(Hp-1) )){
                    // was 10
                    derivCost = 10 * pow( (ZOlGlobal[(i+1)*nz+j] - ZOlGlobal[(i+1+1)*nz+j] ),2); 
                }
                if (( j==2 )&&( i<(Hp-1) )){
                    //was 100
                    derivCost = 100 * pow( (ZOlGlobal[(i+1)*nz+j] - ZOlGlobal[(i+1+1)*nz+j] ),2);
                }
                if (( j==3 )&&( i<(Hp-1) )){
                    derivCost = 0 * pow( (ZOlGlobal[(i+1)*nz+j] - ZOlGlobal[(i+1+1)*nz+j] ),2);
                }
            }

            // ZRefGlobal = given reference trajectory for path following

            /* State cost */
            /* TRANSFERRED *******/
            if  ((trial_number == 1) || ( ZOlGlobal[(i+1)*nz+5]+ s_start - ZRefGlobal[(i+1)*nz+5] > 0 ) ){
                currCostZ = 0.5*Q[j]*pow(ZOlGlobal[(i+1)*nz+j]-ZRefGlobal[(i+1)*nz+j],2);                
                // if lap number = 1 or crossed finish line then cost = path following cost (might be set to 0 in simulation)
            }else{
                if (j == 5){
                    currCostZ = 100;
                }else{
                    currCostZ = 0;
                }
            }
            

            /*if ( (trial_number == 1) || (ReachedSurface == 1) ){
                currCostZ = 0.5*Q[j]*pow(ZOlGlobal[(i+1)*nz+j]-ZRefGlobal[(i+1)*nz+j],2);
            }else if (j == 5){
                    if ( ZOlGlobal[(i+1)*nz+j]+ s_start -ZRefGlobal[(i+1)*nz+j] > 0 ){
                        // if on the surface --> no cost!
                        currCostZ = 0;
                        ReachedSurface = 1;
                    }else{
                        currCostZ = 100;
                    }
            }else{
                    currCostZ = 0;
            }
            */
            
            costZ += currCostZ + derivCost;
            costZSplit[j] += currCostZ;
        }

        /* Lane constraints */
        /* TRANSFERRED ******/
        if (ZOlGlobal[(i+1)*nz+4] > eY_UB)  // soft constraint on upper bound 
        {
            costLane += T[nz-1]*pow((ZOlGlobal[(i+1)*nz+4] - eY_UB), 2);
        }
        else if (ZOlGlobal[(i+1)*nz+4] < eY_LB)  // soft constraint on lower bound 
        {
            costLane += T[nz-1]*pow((eY_LB - ZOlGlobal[(i+1)*nz+4]), 2);
        }
        
    }

    /* Input costs */
    /* TRANSFERRED ******/
    costU = costdU = 0;
    for (i = 0; i < Hu; i++){
        for (j = 0; j < nu; j++){
            costU += 0.5*R[j]*pow(X[i*nu+j] - URefGlobal[i*nu+j],2);
            if (i == 0){
                costdU += 0.5*P[j]*pow(X[i*nu+j]-uPrevGlobal[j],2);     // /pow(dtMPC,2);
            }
            else {
                costdU += 0.5*P[j]*pow(X[i*nu+j]-X[(i-1)*nu+j],2);      // /pow(dt,2);
            }
        }
    }

    /* Slack cost */
    costSlack = 0;
    // if (useQuadSlackCost == 1){
    //     for (i = 0; i < nSlack; i++){
    //         costSlack += T[i]*pow(X[nHu+i],2);
    //     }
    // } else {
    //     for (i = 0; i < nSlack; i++){
    //         costSlack += T[i]*X[nHu+i];
    //     }
    // }

    /* Terminal constraints Initialization */
    /**** TRANSFERRED **********************/
    init_zeros(costZTermCons, nz-1);
    if ( ((trial_number > 1) )){
        for (j = 0; j < nz-2; ++j) //Note this is nz-2 because you don't have to take into account the last state and the filter!
        {
          if ( (trial_number > 2) )
          {                                 // 0 < ParInt < 1 = parameter interpolation
                                            // if we are already in the 3rd lap, we have 2 laps to interpolate from!
                                            // -> we build a new set by the previous interpolations
              termStateErr =       (ParInt*polyval(coeffTermCons[j], nPolyOrderTermCons, ZOlGlobal[Hp*nz+5]) +
                            + (1 - ParInt)*polyval(coeffTermCons_1[j], nPolyOrderTermCons, ZOlGlobal[Hp*nz+5])) - ZOlGlobal[Hp*nz+j];
          }
          else
          {        // in the second lap
              termStateErr =       polyval(coeffTermCons[j], nPolyOrderTermCons, ZOlGlobal[Hp*nz+5]) - ZOlGlobal[Hp*nz+j];  // evaluate polynomial at s
          }
          if (termStateErr > epsTermConsTol[j])     // numerical issue: only add soft constraint cost if it exceeds a certain value (epsTermConsTol)
          {
                costZTermCons[j] += T[j]*pow( (termStateErr - epsTermConsTol[j]) , 2);
          }
          else if (termStateErr < -epsTermConsTol[j])
          {
                costZTermCons[j] +=  T[j]*pow((-epsTermConsTol[j] - termStateErr), 2);
          }
        }
    }
    /* Terminal state cost */
    /**** TRANSFERRED **********************/
    // Q function
    if (ZOlGlobal[Hp*nz+5] + s_start -ZRefGlobal[Hp*nz+5] > 0 ){        // if finish line is crossed
            costZTerm = 0;
    }else{
        if ( (trial_number > 2) ){
             costZTerm =     ParInt*polyval(coeffTermCost, nPolyOrderTermCost, ZOlGlobal[Hp*nz+5]) +
                        + (1-ParInt)*polyval(coeffTermCost_1, nPolyOrderTermCost, ZOlGlobal[Hp*nz+5]);
        }else{
             costZTerm =     polyval(coeffTermCost, nPolyOrderTermCost, ZOlGlobal[Hp*nz+5]);
        }
    }
    
    /* Total cost */
    //*OBJ += costZ + costU + costdU;
    if  ((trial_number == 1) ){ //||( ZOlGlobal[(0+1)*nz+5]+ s_start - ZRefGlobal[(0+1)*nz+5] > 0 ) ){
    //if ( (trial_number == 1) ){
        *OBJ += costZ + costU + costdU + costLane;
    }else{
        *OBJ += costZ + costU + costdU + costSlack + costZTerm +
              + vec_sum(costZTermCons, nz-2) + costLane; //Note this is nz-2 becuase you don't have to take into account the last state and the filter!
    }

    /* Other outputs */
    *MODE = 0;
}

/*** Function implementing the nonlinear constraints for NPSOL ***/
int *funcon(int *MODE, int *ncnln, int *NPARAM, int *NROWJ, int *NEEDC, double *X, double *CON, double *CJAC, int *NSTATE)
{
    int i,j,k;
    double U[nHu];

    /*** Initialization ***/
    for (i = 0; i < nHu; i++){
        U[i] = X[i];
    }

    /*** Simulation ***/
//     simBicycleModel(zCurrGlobal, U);

    /*** Constraints ***/
    /* Dummy */
    CON[0] = 0.0;

    /*** Other outputs ***/
    *MODE = 0;
}

/***** END OF COMPUTATIONAL SOUBROUTINES *****/

/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
	int inpPortWidth[n_Inputs];
    int outPortWidth[n_Outputs];
    int indPort;

	printf("\nINITIALIZING SIZES...");

    /* See sfuntmpl_doc.c for more details on the macros below */

    ssSetNumSFcnParams(S, 0);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    /*** Input ports setup ***/
    if (!ssSetNumInputPorts(S, n_Inputs)) return;
    inpPortWidth[INP_ENABLE]    = 1;
    inpPortWidth[INP_STATE]     = nz;
    inpPortWidth[INP_U_CURR]    = nu;
    inpPortWidth[INP_TUNING]    = nz+2*nu+nSlack+4*nu+6+2;
    inpPortWidth[INP_U_LIN]     = NOPTVAR;
    inpPortWidth[INP_U_REF]     = nu*Hu;        // delta_ref and ax_ref over horizon (accounting for blocking)
    inpPortWidth[INP_Z_REF]     = nz*(Hp+1);    // State reference over horizon
    inpPortWidth[INP_COEFF_TERM_COST] = 2*(nPolyOrderTermCost+1);
    inpPortWidth[INP_COEFF_TERM_CONS] = 2*(nPolyOrderTermCons+1)*(nz-2); // nz-2 beacuse no coeff s and for the filter
    //inpPortWidth[INP_COEFF_TERM_COST] = (nPolyOrderTermCost+1);
    //inpPortWidth[INP_COEFF_TERM_CONS] = (nPolyOrderTermCons+1)*(nz-1);
    inpPortWidth[INP_COEFF_ROAD_ANGLE] = nPolyOrderRoadAngle+1;
    inpPortWidth[INP_COEFF_ROAD_CURV] = nPolyOrderRoadCurv+1;
    inpPortWidth[INP_TRIAL_NUMBER] = 1;
    inpPortWidth[INP_S_START]   = 1;
    inpPortWidth[INP_COEFF_Vy] = 4;
    inpPortWidth[INP_COEFF_YawRate]   = 3;
    inpPortWidth[INP_COEFF_Vx]   = 3;
    
    for (indPort = 0; indPort < n_Inputs; indPort++){
        ssSetInputPortWidth(S, indPort, inpPortWidth[indPort]);
        ssSetInputPortDirectFeedThrough(S, indPort, 1);
        ssSetInputPortRequiredContiguous(S, indPort, 1);
    }

    /*** Output ports setup ***/
    if (!ssSetNumOutputPorts(S, n_Outputs)) return;
	outPortWidth[OUT_DELTA]         = 1;
    outPortWidth[OUT_AX]            = 1;
	outPortWidth[OUT_STAT_OPT]      = 1;
	outPortWidth[OUT_SLACK]         = nSlack;
    outPortWidth[OUT_PARAMS]        = nParamsOut;
    outPortWidth[OUT_U_LIN_S]       = NOPTVAR;
    outPortWidth[OUT_Z_OL]          = nz*(Hp+1);
    outPortWidth[OUT_STATE]         = nz;

    for (indPort = 0; indPort < n_Outputs; indPort++){
        ssSetOutputPortWidth(S, indPort, outPortWidth[indPort]);
    }

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, n_PWorks);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    ssSetOptions(S, 0);
	printf("\n\tSUCCESSFUL!\n");
}


/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    This function is used to specify the sample time(s) for your
 *    S-function. You must register the same number of sample times as
 *    specified in ssSetNumSampleTimes.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, -1);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS   /* Change to #undef to remove function */
#if defined(MDL_INITIALIZE_CONDITIONS)
  /* Function: mdlInitializeConditions ========================================
   * Abstract:
   *    In this function, you should initialize the continuous and discrete
   *    states for your S-function block.  The initial states are placed
   *    in the state vector, ssGetContStates(S) or ssGetRealDiscStates(S).
   *    You can also perform any other initialization activities that your
   *    S-function may require. Note, this routine will be called at the
   *    start of simulation and if it is present in an enabled subsystem
   *    configured to reset states, it will be call when the enabled subsystem
   *    restarts execution to reset the states.
   */
static void mdlInitializeConditions(SimStruct *S)
{

}
#endif /* MDL_INITIALIZE_CONDITIONS */



#define MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START)
  /* Function: mdlStart =======================================================
   * Abstract:
   *    This function is called once at start of model execution. If you
   *    have states that should be initialized once, this is the place
   *    to do it.
   */
static void mdlStart(SimStruct *S)
{

}
#endif /*  MDL_START */



/* Function: mdlOutputs =======================================================
 * Abstract:
 *    In this function, you compute the outputs of your S-function
 *    block. Generally outputs are placed in the output vector, ssGetY(S).
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    /***** DECLARATIONS *****/
    int i, j, k;
    double du_ub[nu], du_lb[nu], u_lb[nu], u_ub[nu];

    /*** Inputs ***/
    int     in_enable    = (int)(((double *) ssGetInputPortSignal(S, INP_ENABLE))[0]);
    double  *in_zCurr    = (double *) ssGetInputPortSignal(S, INP_STATE);
    double  *in_uPrev    = (double *) ssGetInputPortSignal(S, INP_U_CURR);
    double  *in_tuning   = (double *) ssGetInputPortSignal(S, INP_TUNING);
    double  *in_ULin     = (double *) ssGetInputPortSignal(S, INP_U_LIN);
    double  *in_URef     = (double *) ssGetInputPortSignal(S, INP_U_REF);
    double  *in_ZRef     = (double *) ssGetInputPortSignal(S, INP_Z_REF);
    double  *in_coeffTermCost  = (double *) ssGetInputPortSignal(S, INP_COEFF_TERM_COST);
    double  *in_coeffTermCons  = (double *) ssGetInputPortSignal(S, INP_COEFF_TERM_CONS);
    double  *in_coeffRoadAngle = (double *) ssGetInputPortSignal(S, INP_COEFF_ROAD_ANGLE);
    double  *in_coeffRoadCurv  = (double *) ssGetInputPortSignal(S, INP_COEFF_ROAD_CURV);
    int  in_trial_number   = (int)(((double *) ssGetInputPortSignal(S, INP_TRIAL_NUMBER))[0]);
    double  *in_s_start        = (double *) ssGetInputPortSignal(S, INP_S_START);
    double  *in_coeffVy = (double *) ssGetInputPortSignal(S, INP_COEFF_Vy);
    double  *in_coeffYawRate  = (double *) ssGetInputPortSignal(S, INP_COEFF_YawRate);
    double  *in_coeffVx  = (double *) ssGetInputPortSignal(S, INP_COEFF_Vx);
    
    /*** Outputs ***/
	double  *out_delta_f  = (double *) ssGetOutputPortRealSignal(S, OUT_DELTA);
    double  *out_ax       = (double *) ssGetOutputPortRealSignal(S, OUT_AX);
	double  *out_stat_opt = (double *) ssGetOutputPortRealSignal(S, OUT_STAT_OPT);
    double  *out_slack    = (double *) ssGetOutputPortRealSignal(S, OUT_SLACK);
    double  *out_params   = (double *) ssGetOutputPortRealSignal(S, OUT_PARAMS);
    double  *out_ULin     = (double *) ssGetOutputPortRealSignal(S, OUT_U_LIN_S);
    double  *out_ZOl      = (double *) ssGetOutputPortRealSignal(S, OUT_Z_OL);
    double  *out_State    = (double *) ssGetOutputPortRealSignal(S, OUT_STATE);

    /*** NPSOL parameters ***/
    integer noptvar = NOPTVAR;			// number of optimization variables
    integer nclin   = NCLIN;			// number of linear constraints
    integer ncnln   = NCNLN;			// number of non-linear constraints
    integer lda     = NCLIN;			// row dimension of the matrix for linear constraints
    integer ldju    = NCNLN;			// row dimension of the non-linear constraints Jacobian (set to 1 if NCNLN = 0)
    integer ldr     = NCTOTL;			// row dimension of the constraints vector
    integer inform;						// provides info about the solution
    integer iter;						// number of performed iterations
    integer istate[NCTOTL];				// describes the status of the constraints
    doublereal c[NCNLN];				// values of the non-linear constraints
    doublereal cJac[NCNLN][NOPTVAR];	// Jacobian matrix of the non-linear constraints
    doublereal clamda[NCTOTL];			// multipliers of the last solved QP_HL
    doublereal objf = 0.0;				// value of the objective function
    doublereal gradu[NOPTVAR];			// gradient of the objective function
    doublereal Hess[NCTOTL][NOPTVAR];   // Lagrangian Hessian estimate
    doublereal sol[NOPTVAR];			// solution of the problem
    doublereal w[(2*NOPTVAR*NOPTVAR + NOPTVAR*NCLIN + 2*NOPTVAR*NCNLN + 20*NOPTVAR + 11*NCLIN + 21*NCNLN)*2];
    integer lenw = (2*NOPTVAR*NOPTVAR + NOPTVAR*NCLIN + 2*NOPTVAR*NCNLN + 20*NOPTVAR + 11*NCLIN + 21*NCNLN)*2;
    integer iw[(3*NOPTVAR + NCLIN + 2*NCNLN)*2];
    integer leniw = (3*NOPTVAR + NCLIN + 2*NCNLN)*2;

  	/*** Current simulation time ***/
  	t_curr = ssGetT(S);

    /***** GET INPUTS *****/
    /*** Tuning parameters ***/
    /* zCurr, Q */
    init_zeros(Q, nz);
    for (k = 0; k < nz; k++){
        zCurrGlobal[k] = in_zCurr[k];
        Q[k] = in_tuning[k];
    }

    /* R */
    for (k = 0; k < nu; k++){
        R[k] = in_tuning[nz + k];
        P[k] = in_tuning[nz+nu + k];
    }

    /* Slack costs */
    for (k = 0; k < nSlack; k++){
        T[k] = in_tuning[nz+2*nu+k];
    }

    /* Input bounds */
    for (k = 0; k < nu; k++){
        du_lb[k]        = in_tuning[nz+2*nu+nSlack+k];
        du_ub[k]        = in_tuning[nz+3*nu+nSlack+k];
        u_lb[k]         = in_tuning[nz+4*nu+nSlack+k];
        u_ub[k]         = in_tuning[nz+5*nu+nSlack+k];
        uPrevGlobal[k]  = in_uPrev[k];
    }

    /* Miscellaneous tuning parameters */
    useQuadSlackCost = (int) in_tuning[nz+6*nu+nSlack];
    epsTermConsTol[0] = in_tuning[nz+6*nu+nSlack+1];
    epsTermConsTol[1] = in_tuning[nz+6*nu+nSlack+2];
    epsTermConsTol[2] = in_tuning[nz+6*nu+nSlack+3];
    epsTermConsTol[3] = in_tuning[nz+6*nu+nSlack+4];
    epsTermConsTol[4] = in_tuning[nz+6*nu+nSlack+5];

    /* Lane bounds */
    eY_LB = in_tuning[nz+6*nu+nSlack+6];
    eY_UB = in_tuning[nz+6*nu+nSlack+7];

    /*** Driver model inputs ***/
    init_zeros(URefGlobal, nu*Hu);
    for (i = 0; i < Hu; i++){
        for (j = 0; j < nu; j++){
            URefGlobal[i*nu+j] = in_URef[i*nu+j];
        }
    }

   /*** State reference ***/
    for (i = 0; i < Hp+1; i++){
        for (j = 0; j < nz; j++){
            ZRefGlobal[i*nz+j] = in_ZRef[i*nz+j];
        }
    }

    /*** Coefficients ***/
    /* Terminal cost */
    for (i = 0; i <= nPolyOrderTermCost; ++i)
    {
      coeffTermCost[i]   = in_coeffTermCost[i];
      coeffTermCost_1[i] = in_coeffTermCost[(i + nPolyOrderTermCost + 1)];      
    }

    /* Terminal constraint */
    for (j = 0; j < nz-2; ++j) //Note this is nz-2 becuase you don't have to take into account the last state and the filter!
    {
      for (i = 0; i <= nPolyOrderTermCons; ++i)
      {
        coeffTermCons[j][i]   = in_coeffTermCons[j*(nPolyOrderTermCons+1)+i];        
        coeffTermCons_1[j][i] = in_coeffTermCons[(nPolyOrderTermCons+1)*(nz-2) + j*(nPolyOrderTermCons+1)+i];
      }
    }

    /* Road angle */
    for (i = 0; i <= nPolyOrderRoadAngle; ++i)
    {
      coeffRoadAngle[i] = in_coeffRoadAngle[i];
    }

    /* Road curvature */
    for (i = 0; i <= nPolyOrderRoadCurv; ++i)
    {
      coeffRoadCurv[i] = in_coeffRoadCurv[i];
    }

    /* Trial Number*/
    trial_number = in_trial_number;
    s_start = in_s_start[0];

    /* Road curvature */
    for (i = 0; i < 4; ++i)
    {
      coeffVy[i] = in_coeffVy[i];
    }
    /* Road curvature */
    for (i = 0; i < 3; ++i)
    {
      coeffYawRate[i] = in_coeffYawRate[i];
    }
    for (i = 0; i < 3; ++i)
    {
      coeffVx[i] = in_coeffVx[i];
    }
    ReachedSurface = 0;

    /*** Simulation parameters ***/
    init_zeros(out_params,nParamsOut);
    out_params[0] = dt;
    out_params[1] = Hp;
    out_params[2] = Hu;

    /*** Check Enable ***/
	  if (in_enable == 0){
        /* Set default control inputs to driver model inputs */
		    *out_delta_f = 0;     //URefGlobal[0];
		    *out_ax = 0;            //URefGlobal[1];

        *out_stat_opt = -1;
        init_zeros(out_slack,nSlack);

        /* Set default linearization sequence to driver model reference */
        init_zeros(out_ULin, NOPTVAR);

        /* Open loop X,Y */
        simBicycleModel(in_zCurr, out_ULin);
        for (i = 0; i < Hp+1; i++){
            for (j = 0; j < nz; j++){
                out_ZOl[i*nz+j] = ZOlGlobal[i*nz+j];
            }
        }

        /* Other outputs */
        set_const(out_params + 3, nParamsOut - 3, -1);

		    return;
	  }

	/***** INITIALIZATION *****/
    init_zeros(sol,NOPTVAR);
    /*** Initialize solution usfing solution from previous time step ***/
    for (i = 0; i < Hu-1; i++){
        for (j = 0; j < nu; j++){
            sol[i*nu+j] = in_ULin[(i+1)*nu+j];
            if (i == Hu-2){
                sol[(i+1)*nu+j] = in_ULin[(i+1)*nu+j];
            }
        }
    }
    sol[nHu] = in_ULin[nHu];
    /*
     * for (i = 0; i < Hu; i++){
         for (j = 0; j < nu; j++){
             sol[i*nu+j] = in_ULin[(i)*nu+j];
         }
     }
     */
//     
//     /* Initialize ParInt */
//     sol[nHu] = in_ULin[nHu];
    
    /* warm start with previous optimization vector */
//    for (i = 0; i < NOPTVAR; i++){
//        sol[i] = in_ULin[i];
//    }
    
    /***** Solve nonlinear program *****/
    /*** Initialize RHS ***/
    for (k = 0; k < NCTOTL; k++){
        bLGlobal[k] = b_lin_L[k];
        bUGlobal[k] = b_lin_U[k];
    }

    /*** Modify constraint RHS ***/
    /* Input bounds */
    for (i = 0; i < Hu; i++){
        for (j = 0; j < nu; j++){
            bLGlobal[i*nu+j] = u_lb[j];
            bUGlobal[i*nu+j] = u_ub[j];
            if (i==0){
//                 bLGlobal[NOPTVAR+j] = dt*Hi*du_lb[j] + in_uPrev[j];
//                 bUGlobal[NOPTVAR+j] = dt*Hi*du_ub[j] + in_uPrev[j];
                bLGlobal[NOPTVAR+j] = dtMPC*du_lb[j] + in_uPrev[j];
                bUGlobal[NOPTVAR+j] = dtMPC*du_ub[j] + in_uPrev[j];
            }
            else {
                bLGlobal[NOPTVAR+i*nu+j] = dt*Hi*du_lb[j];
                bUGlobal[NOPTVAR+i*nu+j] = dt*Hi*du_ub[j];
            }
        }
    }

    /*** NPSOL subroutines ***/
    npoptn_("Print level = 0", 15);
    npoptn_("Derivative level = 0", 20);
    npoptn_("Optimality Tolerance = 1e-5", 27); //was 5
    npoptn_("Iters = 30", 10); // was 30
    npoptn_("Minor iteration limit = 30", 26);
    npoptn_("Nonlinear feasibility tolerance = 1e-6", 38);

    /*** Call NPSOL ***/
    npsol_(&noptvar, &nclin, &ncnln, &lda, &ldju, &ldr, A_lin, bLGlobal, bUGlobal, funcon, funobj, &inform, &iter,
	     istate, c, cJac, clamda, &objf, gradu, Hess, sol, iw, &leniw, w, &lenw);

    /*** Get optimal values ***/
    *out_delta_f = sol[0];
    *out_ax = sol[1];
//     for (i = 0; i < nSlack; i++){
//         out_slack[i] = sol[nu*Hu+i];
//     }

    /*** Update linearization sequence for next time step ***/
    /* Optimal input sequence */
    for (i = 0; i < Hu; i++){
        for (j = 0; j < nu; j++){
            out_ULin[i*nu+j] = sol[i*nu+j];
        }
    }
    out_ULin[nHu] = sol[nHu];

    /*** Other outputs ***/
    /* Solver outputs */
    *out_stat_opt   = inform;
    out_params[3]   = objf;
    out_params[4]   = iter;
    out_params[5]   = costZ;
    out_params[6]   = costU;
    out_params[7]   = costdU;
    out_params[8]   = costSlack;
    out_params[9]   = vec_max(c,NCNLN);
    for (j = 0; j < nz; j++){
        out_params[10+j] = costZSplit[j];
    }
    out_params[10+nz]  = costLane;
    out_params[10+nz+1]  = costZTerm;
    for (j = 0; j < nz-1; j++){
        out_params[10+nz+2+j]  = costZTermCons[j];
    }
    
    for (j = 0; j < nz-1; j++){
        if (j == 5){
            out_State[j]  = in_zCurr[j] + in_s_start[0];
        }else{
            out_State[j]  = in_zCurr[j];
        }
    }
    
    /* Open loop X, Y */
    simBicycleModel(in_zCurr, sol);
    for (i = 0; i < Hp+1; i++){
        for (j = 0; j < nz; j++){
            out_ZOl[i*nz+j] = ZOlGlobal[i*nz+j];
        }
    }
}

#undef MDL_UPDATE  /* Change to #undef to remove function */
#if defined(MDL_UPDATE)
/* Function: mdlUpdate ======================================================
 * Abstract:
 *    This function is called once for every major integration time step.
   *    Discrete states are typically updated here, but this function is useful
   *    for performing any tasks that should only take place once per
   *    integration step.
   */
static void mdlUpdate(SimStruct *S, int_T tid)
{

}
#endif /* MDL_UPDATE */



#undef MDL_DERIVATIVES  /* Change to #undef to remove function */
#if defined(MDL_DERIVATIVES)
  /* Function: mdlDerivatives =================================================
   * Abstract:
   *    In this function, you compute the S-function block's derivatives.
   *    The derivatives are placed in the derivative vector, ssGetdX(S).
   */
static void mdlDerivatives(SimStruct *S)
{

}
#endif /* MDL_DERIVATIVES */



/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{

}

/*=============================*
 * Required S-function trailer *
 *=============================*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
