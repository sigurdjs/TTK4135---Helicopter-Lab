/*
 * helicopter_no_feedback.c
 *
 * Code generation for model "helicopter_no_feedback".
 *
 * Model version              : 1.207
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Tue Apr 25 10:25:43 2017
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helicopter_no_feedback.h"
#include "helicopter_no_feedback_private.h"
#include "helicopter_no_feedback_dt.h"

/* Block signals (auto storage) */
B_helicopter_no_feedback_T helicopter_no_feedback_B;

/* Continuous states */
X_helicopter_no_feedback_T helicopter_no_feedback_X;

/* Block states (auto storage) */
DW_helicopter_no_feedback_T helicopter_no_feedback_DW;

/* Real-time model */
RT_MODEL_helicopter_no_feedba_T helicopter_no_feedback_M_;
RT_MODEL_helicopter_no_feedba_T *const helicopter_no_feedback_M =
  &helicopter_no_feedback_M_;

/*
 * Writes out MAT-file header.  Returns success or failure.
 * Returns:
 *      0 - success
 *      1 - failure
 */
int_T rt_WriteMat4FileHeader(FILE *fp, int32_T m, int32_T n, const char *name)
{
  typedef enum { ELITTLE_ENDIAN, EBIG_ENDIAN } ByteOrder;

  int16_T one = 1;
  ByteOrder byteOrder = (*((int8_T *)&one)==1) ? ELITTLE_ENDIAN : EBIG_ENDIAN;
  int32_T type = (byteOrder == ELITTLE_ENDIAN) ? 0: 1000;
  int32_T imagf = 0;
  int32_T name_len = (int32_T)strlen(name) + 1;
  if ((fwrite(&type, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&m, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&n, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&imagf, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&name_len, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(name, sizeof(char), name_len, fp) == 0)) {
    return(1);
  } else {
    return(0);
  }
}                                      /* end rt_WriteMat4FileHeader */

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  helicopter_no_feedback_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter_no_feedback_output(void)
{
  /* local block i/o variables */
  real_T rtb_u[2];
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T *lastU;
  real_T rtb_Frontgain;
  real_T rtb_Derivative;
  real_T rtb_Gain1_idx_4;
  real_T rtb_Gain1_idx_5;
  if (rtmIsMajorTimeStep(helicopter_no_feedback_M)) {
    /* set solver stop time */
    if (!(helicopter_no_feedback_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter_no_feedback_M->solverInfo,
                            ((helicopter_no_feedback_M->Timing.clockTickH0 + 1) *
        helicopter_no_feedback_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter_no_feedback_M->solverInfo,
                            ((helicopter_no_feedback_M->Timing.clockTick0 + 1) *
        helicopter_no_feedback_M->Timing.stepSize0 +
        helicopter_no_feedback_M->Timing.clockTickH0 *
        helicopter_no_feedback_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter_no_feedback_M)) {
    helicopter_no_feedback_M->Timing.t[0] = rtsiGetT
      (&helicopter_no_feedback_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helicopter_no_feedback_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S5>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter_no_feedback/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder
        (helicopter_no_feedback_DW.HILReadEncoderTimebase_Task, 1,
         &helicopter_no_feedback_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter_no_feedback_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter_no_feedback_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicopter_no_feedback_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }

    /* Gain: '<S5>/Travel: Count to rad' */
    helicopter_no_feedback_B.TravelCounttorad =
      helicopter_no_feedback_P.TravelCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S12>/Gain' */
    helicopter_no_feedback_B.Gain = helicopter_no_feedback_P.Gain_Gain *
      helicopter_no_feedback_B.TravelCounttorad;

    /* Sum: '<Root>/Sum3' incorporates:
     *  Constant: '<Root>/travel_offset [deg]'
     */
    helicopter_no_feedback_B.Sum3 = helicopter_no_feedback_B.Gain +
      helicopter_no_feedback_P.travel_offsetdeg_Value;
  }

  /* Gain: '<S13>/Gain' incorporates:
   *  TransferFcn: '<S5>/Travel: Transfer Fcn'
   */
  helicopter_no_feedback_B.Gain_d =
    (helicopter_no_feedback_P.TravelTransferFcn_C *
     helicopter_no_feedback_X.TravelTransferFcn_CSTATE +
     helicopter_no_feedback_P.TravelTransferFcn_D *
     helicopter_no_feedback_B.TravelCounttorad) *
    helicopter_no_feedback_P.Gain_Gain_l;
  if (rtmIsMajorTimeStep(helicopter_no_feedback_M)) {
    /* Gain: '<S5>/Pitch: Count to rad' */
    helicopter_no_feedback_B.PitchCounttorad =
      helicopter_no_feedback_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S9>/Gain' */
    helicopter_no_feedback_B.Gain_i = helicopter_no_feedback_P.Gain_Gain_a *
      helicopter_no_feedback_B.PitchCounttorad;
  }

  /* Gain: '<S10>/Gain' incorporates:
   *  TransferFcn: '<S5>/Pitch: Transfer Fcn'
   */
  helicopter_no_feedback_B.Gain_b = (helicopter_no_feedback_P.PitchTransferFcn_C
    * helicopter_no_feedback_X.PitchTransferFcn_CSTATE +
    helicopter_no_feedback_P.PitchTransferFcn_D *
    helicopter_no_feedback_B.PitchCounttorad) *
    helicopter_no_feedback_P.Gain_Gain_ae;
  if (rtmIsMajorTimeStep(helicopter_no_feedback_M)) {
    /* Gain: '<S5>/Elevation: Count to rad' */
    helicopter_no_feedback_B.ElevationCounttorad =
      helicopter_no_feedback_P.ElevationCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S7>/Gain' */
    helicopter_no_feedback_B.Gain_e = helicopter_no_feedback_P.Gain_Gain_lv *
      helicopter_no_feedback_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter_no_feedback_B.Sum = helicopter_no_feedback_B.Gain_e +
      helicopter_no_feedback_P.elavation_offsetdeg_Value;
  }

  /* Gain: '<S8>/Gain' incorporates:
   *  TransferFcn: '<S5>/Elevation: Transfer Fcn'
   */
  helicopter_no_feedback_B.Gain_dg =
    (helicopter_no_feedback_P.ElevationTransferFcn_C *
     helicopter_no_feedback_X.ElevationTransferFcn_CSTATE +
     helicopter_no_feedback_P.ElevationTransferFcn_D *
     helicopter_no_feedback_B.ElevationCounttorad) *
    helicopter_no_feedback_P.Gain_Gain_n;

  /* Gain: '<S3>/Gain1' */
  helicopter_no_feedback_B.Gain1[0] = helicopter_no_feedback_P.Gain1_Gain *
    helicopter_no_feedback_B.Sum3;
  helicopter_no_feedback_B.Gain1[1] = helicopter_no_feedback_P.Gain1_Gain *
    helicopter_no_feedback_B.Gain_d;
  helicopter_no_feedback_B.Gain1[2] = helicopter_no_feedback_P.Gain1_Gain *
    helicopter_no_feedback_B.Gain_i;
  helicopter_no_feedback_B.Gain1[3] = helicopter_no_feedback_P.Gain1_Gain *
    helicopter_no_feedback_B.Gain_b;
  helicopter_no_feedback_B.Gain1[4] = helicopter_no_feedback_P.Gain1_Gain *
    helicopter_no_feedback_B.Sum;
  helicopter_no_feedback_B.Gain1[5] = helicopter_no_feedback_P.Gain1_Gain *
    helicopter_no_feedback_B.Gain_dg;
  if (rtmIsMajorTimeStep(helicopter_no_feedback_M)) {
    /* ToFile: '<Root>/To File' */
    {
      if (!(++helicopter_no_feedback_DW.ToFile_IWORK.Decimation % 1) &&
          (helicopter_no_feedback_DW.ToFile_IWORK.Count*7)+1 < 100000000 ) {
        FILE *fp = (FILE *) helicopter_no_feedback_DW.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[7];
          helicopter_no_feedback_DW.ToFile_IWORK.Decimation = 0;
          u[0] = helicopter_no_feedback_M->Timing.t[1];
          u[1] = helicopter_no_feedback_B.Gain1[0];
          u[2] = helicopter_no_feedback_B.Gain1[1];
          u[3] = helicopter_no_feedback_B.Gain1[2];
          u[4] = helicopter_no_feedback_B.Gain1[3];
          u[5] = helicopter_no_feedback_B.Gain1[4];
          u[6] = helicopter_no_feedback_B.Gain1[5];
          if (fwrite(u, sizeof(real_T), 7, fp) != 7) {
            rtmSetErrorStatus(helicopter_no_feedback_M,
                              "Error writing to MAT-file data.mat");
            return;
          }

          if (((++helicopter_no_feedback_DW.ToFile_IWORK.Count)*7)+1 >=
              100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file data.mat.\n");
          }
        }
      }
    }
  }

  /* FromWorkspace: '<Root>/u+' */
  {
    real_T *pDataValues = (real_T *) helicopter_no_feedback_DW.u_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) helicopter_no_feedback_DW.u_PWORK.TimePtr;
    int_T currTimeIndex = helicopter_no_feedback_DW.u_IWORK.PrevIndex;
    real_T t = helicopter_no_feedback_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[119]) {
      currTimeIndex = 118;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helicopter_no_feedback_DW.u_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_u[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 120;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_u[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 120;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 2; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_u[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 120;
          }
        }
      }
    }
  }

  /* Gain: '<S2>/Gain1' */
  rtb_Gain1_idx_4 = helicopter_no_feedback_P.Gain1_Gain_f *
    helicopter_no_feedback_B.Sum;
  rtb_Gain1_idx_5 = helicopter_no_feedback_P.Gain1_Gain_f *
    helicopter_no_feedback_B.Gain_dg;

  /* Sum: '<S6>/Sum' incorporates:
   *  Constant: '<S6>/Vd_bias'
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S6>/K_pd'
   *  Gain: '<S6>/K_pp'
   *  Sum: '<S6>/Sum2'
   *  Sum: '<S6>/Sum3'
   */
  rtb_Frontgain = ((rtb_u[0] - helicopter_no_feedback_P.Gain1_Gain_f *
                    helicopter_no_feedback_B.Gain_i) *
                   helicopter_no_feedback_P.K_pp -
                   helicopter_no_feedback_P.Gain1_Gain_f *
                   helicopter_no_feedback_B.Gain_b *
                   helicopter_no_feedback_P.K_pd) +
    helicopter_no_feedback_P.Vd_ff;

  /* Integrator: '<S4>/Integrator'
   *
   * Regarding '<S4>/Integrator':
   *  Limited Integrator
   */
  if (helicopter_no_feedback_X.Integrator_CSTATE >=
      helicopter_no_feedback_P.Integrator_UpperSat ) {
    helicopter_no_feedback_X.Integrator_CSTATE =
      helicopter_no_feedback_P.Integrator_UpperSat;
  } else if (helicopter_no_feedback_X.Integrator_CSTATE <=
             (helicopter_no_feedback_P.Integrator_LowerSat) ) {
    helicopter_no_feedback_X.Integrator_CSTATE =
      (helicopter_no_feedback_P.Integrator_LowerSat);
  }

  rtb_Backgain = helicopter_no_feedback_X.Integrator_CSTATE;

  /* Sum: '<S4>/Sum' */
  rtb_Derivative = rtb_u[1] - rtb_Gain1_idx_4;

  /* Sum: '<S4>/Sum2' incorporates:
   *  Constant: '<S4>/Vs_bias'
   *  Gain: '<S4>/K_ed'
   *  Gain: '<S4>/K_ep'
   *  Sum: '<S4>/Sum1'
   */
  rtb_Backgain = ((helicopter_no_feedback_P.K_ep * rtb_Derivative + rtb_Backgain)
                  - helicopter_no_feedback_P.K_ed * rtb_Gain1_idx_5) +
    helicopter_no_feedback_P.Vs_ff;

  /* Sum: '<S1>/Add' */
  rtb_Gain1_idx_4 = rtb_Frontgain + rtb_Backgain;

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (rtb_Backgain - rtb_Frontgain) *
    helicopter_no_feedback_P.Backgain_Gain;

  /* Gain: '<S4>/K_ei' */
  helicopter_no_feedback_B.K_ei = helicopter_no_feedback_P.K_ei * rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter_no_feedback_M)) {
  }

  /* Derivative: '<S5>/Derivative' */
  if ((helicopter_no_feedback_DW.TimeStampA >=
       helicopter_no_feedback_M->Timing.t[0]) &&
      (helicopter_no_feedback_DW.TimeStampB >=
       helicopter_no_feedback_M->Timing.t[0])) {
    rtb_Derivative = 0.0;
  } else {
    rtb_Frontgain = helicopter_no_feedback_DW.TimeStampA;
    lastU = &helicopter_no_feedback_DW.LastUAtTimeA;
    if (helicopter_no_feedback_DW.TimeStampA <
        helicopter_no_feedback_DW.TimeStampB) {
      if (helicopter_no_feedback_DW.TimeStampB <
          helicopter_no_feedback_M->Timing.t[0]) {
        rtb_Frontgain = helicopter_no_feedback_DW.TimeStampB;
        lastU = &helicopter_no_feedback_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter_no_feedback_DW.TimeStampA >=
          helicopter_no_feedback_M->Timing.t[0]) {
        rtb_Frontgain = helicopter_no_feedback_DW.TimeStampB;
        lastU = &helicopter_no_feedback_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helicopter_no_feedback_B.PitchCounttorad - *lastU) /
      (helicopter_no_feedback_M->Timing.t[0] - rtb_Frontgain);
  }

  /* End of Derivative: '<S5>/Derivative' */

  /* Gain: '<S11>/Gain' */
  helicopter_no_feedback_B.Gain_l = helicopter_no_feedback_P.Gain_Gain_a1 *
    rtb_Derivative;
  if (rtmIsMajorTimeStep(helicopter_no_feedback_M)) {
  }

  /* Saturate: '<S5>/Back motor: Saturation' */
  if (rtb_Backgain > helicopter_no_feedback_P.BackmotorSaturation_UpperSat) {
    helicopter_no_feedback_B.BackmotorSaturation =
      helicopter_no_feedback_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain <
             helicopter_no_feedback_P.BackmotorSaturation_LowerSat) {
    helicopter_no_feedback_B.BackmotorSaturation =
      helicopter_no_feedback_P.BackmotorSaturation_LowerSat;
  } else {
    helicopter_no_feedback_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S5>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_no_feedback_M)) {
  }

  /* Gain: '<S1>/Front gain' */
  rtb_Gain1_idx_4 *= helicopter_no_feedback_P.Frontgain_Gain;

  /* Saturate: '<S5>/Front motor: Saturation' */
  if (rtb_Gain1_idx_4 > helicopter_no_feedback_P.FrontmotorSaturation_UpperSat)
  {
    helicopter_no_feedback_B.FrontmotorSaturation =
      helicopter_no_feedback_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Gain1_idx_4 <
             helicopter_no_feedback_P.FrontmotorSaturation_LowerSat) {
    helicopter_no_feedback_B.FrontmotorSaturation =
      helicopter_no_feedback_P.FrontmotorSaturation_LowerSat;
  } else {
    helicopter_no_feedback_B.FrontmotorSaturation = rtb_Gain1_idx_4;
  }

  /* End of Saturate: '<S5>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_no_feedback_M)) {
    /* S-Function (hil_write_analog_block): '<S5>/HIL Write Analog' */

    /* S-Function Block: helicopter_no_feedback/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter_no_feedback_DW.HILWriteAnalog_Buffer[0] =
        helicopter_no_feedback_B.FrontmotorSaturation;
      helicopter_no_feedback_DW.HILWriteAnalog_Buffer[1] =
        helicopter_no_feedback_B.BackmotorSaturation;
      result = hil_write_analog(helicopter_no_feedback_DW.HILInitialize_Card,
        helicopter_no_feedback_P.HILWriteAnalog_channels, 2,
        &helicopter_no_feedback_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helicopter_no_feedback_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S5>/Derivative' */
  if (helicopter_no_feedback_DW.TimeStampA == (rtInf)) {
    helicopter_no_feedback_DW.TimeStampA = helicopter_no_feedback_M->Timing.t[0];
    lastU = &helicopter_no_feedback_DW.LastUAtTimeA;
  } else if (helicopter_no_feedback_DW.TimeStampB == (rtInf)) {
    helicopter_no_feedback_DW.TimeStampB = helicopter_no_feedback_M->Timing.t[0];
    lastU = &helicopter_no_feedback_DW.LastUAtTimeB;
  } else if (helicopter_no_feedback_DW.TimeStampA <
             helicopter_no_feedback_DW.TimeStampB) {
    helicopter_no_feedback_DW.TimeStampA = helicopter_no_feedback_M->Timing.t[0];
    lastU = &helicopter_no_feedback_DW.LastUAtTimeA;
  } else {
    helicopter_no_feedback_DW.TimeStampB = helicopter_no_feedback_M->Timing.t[0];
    lastU = &helicopter_no_feedback_DW.LastUAtTimeB;
  }

  *lastU = helicopter_no_feedback_B.PitchCounttorad;

  /* End of Update for Derivative: '<S5>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter_no_feedback_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter_no_feedback_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helicopter_no_feedback_M->Timing.clockTick0)) {
    ++helicopter_no_feedback_M->Timing.clockTickH0;
  }

  helicopter_no_feedback_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helicopter_no_feedback_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.002s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helicopter_no_feedback_M->Timing.clockTick1)) {
      ++helicopter_no_feedback_M->Timing.clockTickH1;
    }

    helicopter_no_feedback_M->Timing.t[1] =
      helicopter_no_feedback_M->Timing.clockTick1 *
      helicopter_no_feedback_M->Timing.stepSize1 +
      helicopter_no_feedback_M->Timing.clockTickH1 *
      helicopter_no_feedback_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter_no_feedback_derivatives(void)
{
  XDot_helicopter_no_feedback_T *_rtXdot;
  _rtXdot = ((XDot_helicopter_no_feedback_T *)
             helicopter_no_feedback_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S5>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE +=
    helicopter_no_feedback_P.TravelTransferFcn_A *
    helicopter_no_feedback_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter_no_feedback_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S5>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE +=
    helicopter_no_feedback_P.PitchTransferFcn_A *
    helicopter_no_feedback_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter_no_feedback_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S5>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter_no_feedback_P.ElevationTransferFcn_A *
    helicopter_no_feedback_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE +=
    helicopter_no_feedback_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S4>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helicopter_no_feedback_X.Integrator_CSTATE <=
            (helicopter_no_feedback_P.Integrator_LowerSat) );
    usat = ( helicopter_no_feedback_X.Integrator_CSTATE >=
            helicopter_no_feedback_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helicopter_no_feedback_B.K_ei > 0)) ||
        (usat && (helicopter_no_feedback_B.K_ei < 0)) ) {
      ((XDot_helicopter_no_feedback_T *)
        helicopter_no_feedback_M->ModelData.derivs)->Integrator_CSTATE =
        helicopter_no_feedback_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helicopter_no_feedback_T *)
        helicopter_no_feedback_M->ModelData.derivs)->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helicopter_no_feedback_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter_no_feedback/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0",
                      &helicopter_no_feedback_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options
      (helicopter_no_feedback_DW.HILInitialize_Card,
       "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter_no_feedback_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
      return;
    }

    if ((helicopter_no_feedback_P.HILInitialize_set_analog_input_ &&
         !is_switching) ||
        (helicopter_no_feedback_P.HILInitialize_set_analog_inpu_m &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums =
          &helicopter_no_feedback_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] =
            helicopter_no_feedback_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums =
          &helicopter_no_feedback_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] =
            helicopter_no_feedback_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges
        (helicopter_no_feedback_DW.HILInitialize_Card,
         helicopter_no_feedback_P.HILInitialize_analog_input_chan, 8U,
         &helicopter_no_feedback_DW.HILInitialize_AIMinimums[0],
         &helicopter_no_feedback_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_no_feedback_P.HILInitialize_set_analog_output &&
         !is_switching) ||
        (helicopter_no_feedback_P.HILInitialize_set_analog_outp_b &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums =
          &helicopter_no_feedback_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] =
            helicopter_no_feedback_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums =
          &helicopter_no_feedback_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] =
            helicopter_no_feedback_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges
        (helicopter_no_feedback_DW.HILInitialize_Card,
         helicopter_no_feedback_P.HILInitialize_analog_output_cha, 8U,
         &helicopter_no_feedback_DW.HILInitialize_AOMinimums[0],
         &helicopter_no_feedback_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_no_feedback_P.HILInitialize_set_analog_outp_e &&
         !is_switching) ||
        (helicopter_no_feedback_P.HILInitialize_set_analog_outp_j &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages =
          &helicopter_no_feedback_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter_no_feedback_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helicopter_no_feedback_DW.HILInitialize_Card,
        helicopter_no_feedback_P.HILInitialize_analog_output_cha, 8U,
        &helicopter_no_feedback_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_no_feedback_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages =
          &helicopter_no_feedback_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter_no_feedback_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter_no_feedback_DW.HILInitialize_Card,
         helicopter_no_feedback_P.HILInitialize_analog_output_cha, 8U,
         &helicopter_no_feedback_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_no_feedback_P.HILInitialize_set_encoder_param &&
         !is_switching) ||
        (helicopter_no_feedback_P.HILInitialize_set_encoder_par_m &&
         is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter_no_feedback_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] =
            helicopter_no_feedback_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helicopter_no_feedback_DW.HILInitialize_Card,
         helicopter_no_feedback_P.HILInitialize_encoder_channels, 8U,
         (t_encoder_quadrature_mode *)
         &helicopter_no_feedback_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_no_feedback_P.HILInitialize_set_encoder_count &&
         !is_switching) ||
        (helicopter_no_feedback_P.HILInitialize_set_encoder_cou_k &&
         is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter_no_feedback_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] =
            helicopter_no_feedback_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts
        (helicopter_no_feedback_DW.HILInitialize_Card,
         helicopter_no_feedback_P.HILInitialize_encoder_channels, 8U,
         &helicopter_no_feedback_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_no_feedback_P.HILInitialize_set_pwm_params_at &&
         !is_switching) ||
        (helicopter_no_feedback_P.HILInitialize_set_pwm_params__f &&
         is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter_no_feedback_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter_no_feedback_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helicopter_no_feedback_DW.HILInitialize_Card,
        helicopter_no_feedback_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helicopter_no_feedback_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helicopter_no_feedback_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues =
          &helicopter_no_feedback_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helicopter_no_feedback_DW.HILInitialize_POSortedChans[num_duty_cycle_modes]
              = p_HILInitialize_pwm_channels[i1];
            helicopter_no_feedback_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]
              = helicopter_no_feedback_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helicopter_no_feedback_DW.HILInitialize_POSortedChans[7U -
              num_frequency_modes] = p_HILInitialize_pwm_channels[i1];
            helicopter_no_feedback_DW.HILInitialize_POSortedFreqs[7U -
              num_frequency_modes] =
              helicopter_no_feedback_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency
          (helicopter_no_feedback_DW.HILInitialize_Card,
           &helicopter_no_feedback_DW.HILInitialize_POSortedChans[0],
           num_duty_cycle_modes,
           &helicopter_no_feedback_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle
          (helicopter_no_feedback_DW.HILInitialize_Card,
           &helicopter_no_feedback_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
           num_frequency_modes,
           &helicopter_no_feedback_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues =
          &helicopter_no_feedback_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] =
            helicopter_no_feedback_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helicopter_no_feedback_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] =
            helicopter_no_feedback_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter_no_feedback_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] =
            helicopter_no_feedback_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration
        (helicopter_no_feedback_DW.HILInitialize_Card,
         helicopter_no_feedback_P.HILInitialize_pwm_channels, 8U,
         (t_pwm_configuration *)
         &helicopter_no_feedback_DW.HILInitialize_POModeValues[0],
         (t_pwm_alignment *)
         &helicopter_no_feedback_DW.HILInitialize_POAlignValues[0],
         (t_pwm_polarity *)
         &helicopter_no_feedback_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs =
          &helicopter_no_feedback_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] =
            helicopter_no_feedback_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter_no_feedback_DW.HILInitialize_POValues
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] =
            helicopter_no_feedback_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helicopter_no_feedback_DW.HILInitialize_Card,
        helicopter_no_feedback_P.HILInitialize_pwm_channels, 8U,
        &helicopter_no_feedback_DW.HILInitialize_POSortedFreqs[0],
        &helicopter_no_feedback_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_no_feedback_P.HILInitialize_set_pwm_outputs_a &&
         !is_switching) ||
        (helicopter_no_feedback_P.HILInitialize_set_pwm_outputs_g &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_no_feedback_DW.HILInitialize_POValues
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] =
            helicopter_no_feedback_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helicopter_no_feedback_DW.HILInitialize_Card,
        helicopter_no_feedback_P.HILInitialize_pwm_channels, 8U,
        &helicopter_no_feedback_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_no_feedback_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_no_feedback_DW.HILInitialize_POValues
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] =
            helicopter_no_feedback_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter_no_feedback_DW.HILInitialize_Card,
         helicopter_no_feedback_P.HILInitialize_pwm_channels, 8U,
         &helicopter_no_feedback_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S5>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter_no_feedback/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (helicopter_no_feedback_DW.HILInitialize_Card,
       helicopter_no_feedback_P.HILReadEncoderTimebase_samples_,
       helicopter_no_feedback_P.HILReadEncoderTimebase_channels, 3,
       &helicopter_no_feedback_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
    }
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    char fileName[509] = "data.mat";
    FILE *fp = (NULL);
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helicopter_no_feedback_M,
                        "Error creating .mat file data.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,7,0,"data")) {
      rtmSetErrorStatus(helicopter_no_feedback_M,
                        "Error writing mat file header to file data.mat");
      return;
    }

    helicopter_no_feedback_DW.ToFile_IWORK.Count = 0;
    helicopter_no_feedback_DW.ToFile_IWORK.Decimation = -1;
    helicopter_no_feedback_DW.ToFile_PWORK.FilePtr = fp;
  }

  /* Start for FromWorkspace: '<Root>/u+' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.49481903599531779, 0.16014638821386434,
      -0.10626328803862102, -0.307278306480395, -0.45154407399349744,
      -0.52359877559829882, -0.52359877559829882, -0.52359877559829882,
      -0.52359877559829882, -0.52359877559829882, -0.52359877559829882,
      -0.51159896753470135, -0.45712995638768061, -0.39980093911198422,
      -0.34221206095736473, -0.28636565715661072, -0.23375737331410953,
      -0.18542470303039477, -0.14201485766965574, -0.10385960931000998,
      -0.071039027351489342, -0.043427010736115788, -0.020760545532702145,
      -0.0026532057262954349, 0.011331646956018356, 0.021677856286943045,
      0.028879262711539077, 0.033421003337928468, 0.035766054664039808,
      0.036347684240999212, 0.035555194730558237, 0.033731021427047821,
      0.031183024678530902, 0.028164528270653449, 0.024885650622257477,
      0.02152182668138166, 0.018205129131723567, 0.015040281409484685,
      0.012097110072578318, 0.0094268696131999154, 0.0070610468302674465,
      0.0049992382690393113, 0.0032498204497654875, 0.0017986264491737466,
      0.00062615902922036142, -0.00028966582250303751, -0.00097947374119644061,
      -0.0014745362668541287, -0.0017935498295042751, -0.0019725952606144985,
      -0.0020364085801468051, -0.0020027877995300422, -0.0019027113680988479,
      -0.0017540636979914007, -0.0015693552767636067, -0.0013615937826847641,
      -0.0011493378349797083, -0.00093921374088163192, -0.00073972138758464649,
      -0.00055603764814935331, -0.00039817301995658, -0.00026478930263344963,
      -0.00016007725404975484, -8.39090499064568E-5, -3.5329005868231011E-5,
      -8.9887773023080052E-6, 1.3919686989920261E-6, 1.4510386305556906E-6,
      3.6870306124610452E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.71545654317320562, 0.75937650943698232, 0.80103474243544281,
      0.8391501086969172, 0.87183571027071494, 0.89675907433319835,
      0.91101279639517951, 0.91102362520874258, 0.89246904997035581,
      0.84996086965642537, 0.77711762179486477, 0.66611874081333888,
      0.507495313990935, 0.29001235159535238, 2.8737646339860942E-5,
      1.6370287059902176E-5, -1.0143382896995259E-5, -2.6971281605155993E-5,
      -2.1695338889045957E-5, -4.3996246426019175E-5, -3.4924126271687144E-5,
      -2.1589015641204992E-5, -3.2157272187417113E-5, -2.1992910839530436E-5,
      -2.7363799948772867E-5, -6.6953186066275473E-6, -1.1898235305523046E-5,
      -1.5404046412305385E-6, 3.5232933716780842E-6, 1.5503829315593544E-5,
      1.357637826462154E-5, 6.4322730594583778E-6, 1.9422004002202586E-6,
      3.3293660113314164E-5, 2.5189956546357963E-5, 3.8724705811281176E-5,
      2.8878215630788409E-5, 3.4199995271265545E-5, 2.8093635584241423E-5,
      1.6165759695328223E-5, 3.2679416447857613E-5, 2.093622778025257E-5,
      1.4858823365536103E-5, 3.1667354646397814E-5, 2.6237545660898295E-5,
      1.3198010435862133E-5, 1.68712302557705E-5, 2.0710787406506706E-5,
      1.1034642451954818E-5, 1.4315016627632981E-5, 6.2594324714316339E-6,
      8.327067860064766E-6, 7.363440181698258E-6, 5.3953379464852881E-6,
      -1.1184252033667726E-5, -1.4676824255462339E-5, -9.7597820396959535E-6,
      -2.1649270797334584E-6, -3.400609951150895E-6, 7.5845213698175247E-6,
      6.7875543635065915E-6, -1.0563353219173227E-5, -3.7022187473625117E-6,
      1.88209906733747E-6, -1.3809070078934608E-5, -1.02233183817887E-5,
      -5.4570463715209576E-6, 3.6151058685247505E-7, 3.5835016919799951E-6,
      3.799435053319569E-6, -1.2240198287727813E-5, 3.0392581368119163E-7,
      1.0424528015192125E-6, 4.7439105154721909E-6, -8.6495195628449869E-6,
      -4.1145733636697511E-6, 1.7617853831570765E-5, 1.7925046229036758E-5,
      -6.0065290702303496E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter_no_feedback_DW.u_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter_no_feedback_DW.u_PWORK.DataPtr = (void *) pDataValues0;
    helicopter_no_feedback_DW.u_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S5>/Travel: Transfer Fcn' */
  helicopter_no_feedback_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S5>/Pitch: Transfer Fcn' */
  helicopter_no_feedback_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S5>/Elevation: Transfer Fcn' */
  helicopter_no_feedback_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S4>/Integrator' */
  helicopter_no_feedback_X.Integrator_CSTATE =
    helicopter_no_feedback_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S5>/Derivative' */
  helicopter_no_feedback_DW.TimeStampA = (rtInf);
  helicopter_no_feedback_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter_no_feedback_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter_no_feedback/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter_no_feedback_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter_no_feedback_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter_no_feedback_P.HILInitialize_set_analog_out_ex &&
         !is_switching) ||
        (helicopter_no_feedback_P.HILInitialize_set_analog_outp_c &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages =
          &helicopter_no_feedback_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] =
            helicopter_no_feedback_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helicopter_no_feedback_P.HILInitialize_set_pwm_output_ap &&
         !is_switching) ||
        (helicopter_no_feedback_P.HILInitialize_set_pwm_outputs_p &&
         is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_no_feedback_DW.HILInitialize_POValues
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] =
            helicopter_no_feedback_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter_no_feedback_DW.HILInitialize_Card
                         ,
                         helicopter_no_feedback_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helicopter_no_feedback_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter_no_feedback_DW.HILInitialize_AOVoltages[0]
                         , &helicopter_no_feedback_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog
            (helicopter_no_feedback_DW.HILInitialize_Card,
             helicopter_no_feedback_P.HILInitialize_analog_output_cha,
             num_final_analog_outputs,
             &helicopter_no_feedback_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm
            (helicopter_no_feedback_DW.HILInitialize_Card,
             helicopter_no_feedback_P.HILInitialize_pwm_channels,
             num_final_pwm_outputs,
             &helicopter_no_feedback_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_no_feedback_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter_no_feedback_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter_no_feedback_DW.HILInitialize_Card);
    hil_close(helicopter_no_feedback_DW.HILInitialize_Card);
    helicopter_no_feedback_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helicopter_no_feedback_DW.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "data.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_no_feedback_M,
                          "Error closing MAT-file data.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helicopter_no_feedback_M,
                          "Error reopening MAT-file data.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 7,
           helicopter_no_feedback_DW.ToFile_IWORK.Count, "data")) {
        rtmSetErrorStatus(helicopter_no_feedback_M,
                          "Error writing header for data to MAT-file data.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_no_feedback_M,
                          "Error closing MAT-file data.mat");
        return;
      }

      helicopter_no_feedback_DW.ToFile_PWORK.FilePtr = (NULL);
    }
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  helicopter_no_feedback_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter_no_feedback_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  helicopter_no_feedback_initialize();
}

void MdlTerminate(void)
{
  helicopter_no_feedback_terminate();
}

/* Registration function */
RT_MODEL_helicopter_no_feedba_T *helicopter_no_feedback(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter_no_feedback_P.Integrator_UpperSat = rtInf;
  helicopter_no_feedback_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter_no_feedback_M, 0,
                sizeof(RT_MODEL_helicopter_no_feedba_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter_no_feedback_M->solverInfo,
                          &helicopter_no_feedback_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter_no_feedback_M->solverInfo, &rtmGetTPtr
                (helicopter_no_feedback_M));
    rtsiSetStepSizePtr(&helicopter_no_feedback_M->solverInfo,
                       &helicopter_no_feedback_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter_no_feedback_M->solverInfo,
                 &helicopter_no_feedback_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helicopter_no_feedback_M->solverInfo, (real_T **)
                         &helicopter_no_feedback_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helicopter_no_feedback_M->solverInfo,
      &helicopter_no_feedback_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helicopter_no_feedback_M->solverInfo,
                          (&rtmGetErrorStatus(helicopter_no_feedback_M)));
    rtsiSetRTModelPtr(&helicopter_no_feedback_M->solverInfo,
                      helicopter_no_feedback_M);
  }

  rtsiSetSimTimeStep(&helicopter_no_feedback_M->solverInfo, MAJOR_TIME_STEP);
  helicopter_no_feedback_M->ModelData.intgData.f[0] =
    helicopter_no_feedback_M->ModelData.odeF[0];
  helicopter_no_feedback_M->ModelData.contStates = ((real_T *)
    &helicopter_no_feedback_X);
  rtsiSetSolverData(&helicopter_no_feedback_M->solverInfo, (void *)
                    &helicopter_no_feedback_M->ModelData.intgData);
  rtsiSetSolverName(&helicopter_no_feedback_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter_no_feedback_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter_no_feedback_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter_no_feedback_M->Timing.sampleTimes =
      (&helicopter_no_feedback_M->Timing.sampleTimesArray[0]);
    helicopter_no_feedback_M->Timing.offsetTimes =
      (&helicopter_no_feedback_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter_no_feedback_M->Timing.sampleTimes[0] = (0.0);
    helicopter_no_feedback_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helicopter_no_feedback_M->Timing.offsetTimes[0] = (0.0);
    helicopter_no_feedback_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter_no_feedback_M, &helicopter_no_feedback_M->Timing.tArray
             [0]);

  {
    int_T *mdlSampleHits = helicopter_no_feedback_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter_no_feedback_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter_no_feedback_M, -1);
  helicopter_no_feedback_M->Timing.stepSize0 = 0.002;
  helicopter_no_feedback_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helicopter_no_feedback_M->Sizes.checksums[0] = (1132350995U);
  helicopter_no_feedback_M->Sizes.checksums[1] = (3543923485U);
  helicopter_no_feedback_M->Sizes.checksums[2] = (3882204645U);
  helicopter_no_feedback_M->Sizes.checksums[3] = (1401270199U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helicopter_no_feedback_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helicopter_no_feedback_M->extModeInfo,
      &helicopter_no_feedback_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter_no_feedback_M->extModeInfo,
                        helicopter_no_feedback_M->Sizes.checksums);
    rteiSetTPtr(helicopter_no_feedback_M->extModeInfo, rtmGetTPtr
                (helicopter_no_feedback_M));
  }

  helicopter_no_feedback_M->solverInfoPtr =
    (&helicopter_no_feedback_M->solverInfo);
  helicopter_no_feedback_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter_no_feedback_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter_no_feedback_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter_no_feedback_M->ModelData.blockIO = ((void *)
    &helicopter_no_feedback_B);

  {
    int_T i;
    for (i = 0; i < 6; i++) {
      helicopter_no_feedback_B.Gain1[i] = 0.0;
    }

    helicopter_no_feedback_B.TravelCounttorad = 0.0;
    helicopter_no_feedback_B.Gain = 0.0;
    helicopter_no_feedback_B.Sum3 = 0.0;
    helicopter_no_feedback_B.Gain_d = 0.0;
    helicopter_no_feedback_B.PitchCounttorad = 0.0;
    helicopter_no_feedback_B.Gain_i = 0.0;
    helicopter_no_feedback_B.Gain_b = 0.0;
    helicopter_no_feedback_B.ElevationCounttorad = 0.0;
    helicopter_no_feedback_B.Gain_e = 0.0;
    helicopter_no_feedback_B.Sum = 0.0;
    helicopter_no_feedback_B.Gain_dg = 0.0;
    helicopter_no_feedback_B.K_ei = 0.0;
    helicopter_no_feedback_B.Gain_l = 0.0;
    helicopter_no_feedback_B.BackmotorSaturation = 0.0;
    helicopter_no_feedback_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helicopter_no_feedback_M->ModelData.defaultParam = ((real_T *)
    &helicopter_no_feedback_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter_no_feedback_X;
    helicopter_no_feedback_M->ModelData.contStates = (x);
    (void) memset((void *)&helicopter_no_feedback_X, 0,
                  sizeof(X_helicopter_no_feedback_T));
  }

  /* states (dwork) */
  helicopter_no_feedback_M->ModelData.dwork = ((void *)
    &helicopter_no_feedback_DW);
  (void) memset((void *)&helicopter_no_feedback_DW, 0,
                sizeof(DW_helicopter_no_feedback_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_no_feedback_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_no_feedback_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_no_feedback_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_no_feedback_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_no_feedback_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_no_feedback_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_no_feedback_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helicopter_no_feedback_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter_no_feedback_DW.TimeStampA = 0.0;
  helicopter_no_feedback_DW.LastUAtTimeA = 0.0;
  helicopter_no_feedback_DW.TimeStampB = 0.0;
  helicopter_no_feedback_DW.LastUAtTimeB = 0.0;
  helicopter_no_feedback_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter_no_feedback_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter_no_feedback_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter_no_feedback_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter_no_feedback_M->Sizes.numY = (0);/* Number of model outputs */
  helicopter_no_feedback_M->Sizes.numU = (0);/* Number of model inputs */
  helicopter_no_feedback_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter_no_feedback_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter_no_feedback_M->Sizes.numBlocks = (54);/* Number of blocks */
  helicopter_no_feedback_M->Sizes.numBlockIO = (16);/* Number of block outputs */
  helicopter_no_feedback_M->Sizes.numBlockPrms = (142);/* Sum of parameter "widths" */
  return helicopter_no_feedback_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
