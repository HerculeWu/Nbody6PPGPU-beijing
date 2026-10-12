*     TidalTiles pc-Msun-Myr-v1: match PeTar -u 1 and frozen analysis.
*     One source for SCALE, UNITS and MWPotential2014. Not a conversion
*     for legacy snapshots/checkpoints: those retain their old scales.
      REAL*8 TT_KMS_PCMYR,TT_G_PCMYR,TT_G_KMS,TT_PC_CM
      PARAMETER (TT_KMS_PCMYR=1.022712165045695D0)
      PARAMETER (TT_G_PCMYR=0.00449830997959438D0)
      PARAMETER (TT_G_KMS=TT_G_PCMYR/TT_KMS_PCMYR**2)
*     Julian Myr = 31557600 * 10^6 seconds; one km = 10^5 cm.
      PARAMETER (TT_PC_CM=3.15576D13*1.0D5/TT_KMS_PCMYR)
      CHARACTER*14 TT_UNIT_VERSION
      PARAMETER (TT_UNIT_VERSION='pc-msun-myr-v1')
