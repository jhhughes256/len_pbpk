$PROB LENA PBPK
$INPUT ID TIME AMT EVID DV CMT MDV WT TADNOM DOSEMGKG DOSE ROUTE
$DATA ..\nmprep.csv IGNORE=C

$SUBR ADVAN8 TOL=2 ; DES stiff

$MODEL
  COMP = (PA)  ; 1 vascular mixing
  COMP = (ART)  ; 2 arterial blood
  COMP = (BRA)  ; 3 brain
  COMP = (LVR)  ; 4 liver
  COMP = (GUT)  ; 5 gut
  COMP = (SPLR)  ; 6 spleen sinus
  COMP = (SPLS)  ; 7 spleen pulp
  COMP = (KID)  ; 8 kidney
  COMP = (FILT)  ; 9 renal tubule filtrate
  COMP = (TUBC)  ; 10 renal tubule cells
  COMP = (URINE)  ; 11 urine
  COMP = (HRT)  ; 12 heart
  COMP = (MSC)  ; 13 muscle
  COMP = (BOD)  ; 14 rest of body
  COMP = (INT)  ; 15 oral depot

$PK
; standard physiological parameters
  WT_STD=25 ; weight (g)
  CO_STD=13.98 ; cardiac output (ml/min)
  
; cardiac output and volume for rest of body
  QBOD_STD=CO_STD-(QBRA_STD+QLVR_STD+QKID_STD+QHRT_STD+QMSC_STD)
  VBOD_STD=WT_STD-(VBRA_STD+VLVR_STD+VGUT_STD+VSPL_STD+VKID_STD+VHRT_STD+VMSC_STD)

; allometric scaling of blood flows
  QBRA=QBRA_STD*(WT/WT_STD)**0.75
  QLVR=QLVR_STD*(WT/WT_STD)**0.75
  QGUT=QGUT_STD*(WT/WT_STD)**0.75
  QSPL=QSPL_STD*(WT/WT_STD)**0.75
  QKID=QKID_STD*(WT/WT_STD)**0.75
  QHRT=QHRT_STD*(WT/WT_STD)**0.75
  QMSC=QMSC_STD*(WT/WT_STD)**0.75
  QBOD=QBOD_STD*(WT/WT_STD)**0.75
  QCO=QBRA+QLVR+QKID+QHRT+QMSC+QBOD
  
; allometric scaling for clearances and permeabilities
  PSSPL=PSSPL_STD*(WT/WT_STD)**0.75
  GFR=GFR_STD*(WT/WT_STD)**0.75
  KURINE=KURINE_STD*(WT/WT_STD)**0.75
  PSDIFF=PSDIFF_STD*(WT/WT_STD)**0.75
  KA=KA_STD
  
; apparent distribution volumes with allometric scaling
  VMIX=VMIX_STD*EXP(PPVFUT)*(WT/WT_STD)**1  ; no apparent distribution component
  VLNG=VLNG_STD*(FU/FUTLNG)*(WT/WT_STD)**1
  VBRA=VBRA_STD*(FU/FUTBRA)*(WT/WT_STD)**1
  VLVR=VLVR_STD*(FU/FUTLVR)*(WT/WT_STD)**1
  VGUT=VGUT_STD*(FU/FUTGUT)*(WT/WT_STD)**1
  VSPL=VSPL_STD*(FU/FUTSPL)*(WT/WT_STD)**1
  VKID=VKID_STD*(FU/FUTKID)*(WT/WT_STD)**1
  VHRT=VHRT_STD*(FU/FUTHRT)*(WT/WT_STD)**1
  VMSC=VMSC_STD*(FU/FUTMSC)*(WT/WT_STD)**1
  VBOD=VBOD_STD*(FU/FUTBOD)*(WT/WT_STD)**1

$DES
  DADT(1) = QBRA*A(3)/VBRA +QLVR*A(4)/VLVR +QKID*A(8)/VKID +QHRT*A(12)/VHRT +QMSC*A(13)/VMSC +QBOD*A(14)/VBOD -QCO*A(1)/VMIX  ; vascular mixing
  DADT(2) = QCO*(A(1)/VMIX -A(2)/VLNG)  ; arterial blood
  DADT(3) = QBRA*(A(2)/VLNG -A(3)/VBRA)  ; brain
  DADT(4) = (QLVR-(QSPL+QGUT))*A(2)/VLNG +QGUT*A(5)/VGUT +QSPL*A(6)/VSPL -QLVR*A(4)/VLVR  ; liver
  DADT(5) = QGUT*(A(2)/VLNG -A(5)/VGUT) +KA*A(15)  ; gut
  DADT(6) = QSPL*(A(2)/VLNG -A(6)/VSPL) +PSSPL*(A(7) -A(6))  ; spleen sinus
  DADT(7) = PSSPL*(A(6) -A(7))  ; spleen pulp
  DADT(8) = QKID*(A(2)/VLNG -A(8)/VKID) -GFR*A(2)/VLNG +PSDIFF*(A(10) -A(8))  ; kidney
  PSTRAN = VMAX/(KM + A(8)/VKID)  ; saturable excretion
  DADT(9) = GFR*A(2)/VLNG +PSTRAN*A(10) -KURINE*A(9)  ; renal tubule filtrate
  DADT(10) = PSDIFF*(A(8) -A(10)) -PSTRAN*A(10)  ; renal tubule cells
  DADT(11) = KURINE*A(9)  ; urine
  DADT(12) = QHRT*(A(2)/VLNG -A(12)/VHRT)  ; heart
  DADT(13) = QMSC*(A(2)/VLNG -A(13)/VMSC)  ; muscle
  DADT(14) = QBOD*(A(2)/VLNG -A(14)/VBOD)  ; rest of body
  DADT(15) = -KA*A(15)

$ERROR
  A_PA = A(1)
  A_LNG = A(2)
  A_BRA = A(3)
  A_LVR = A(4)
  A_GUT = A(5)
  A_SPLR = A(6)
  A_SPLS = A(7)
  A_KID = A(8)
  A_FILT = A(9)
  A_TUBC = A(10)
  A_URINE = A(11)
  A_HRT = A(12)
  A_MSC = A(13)
  A_BOD = A(14)
  A_PO = A(15)

  PRED_PA = A_PA/VMIX
  PRED_LNG = A_LNG/VLNG
  PRED_BRA = A_BRA/VBRA
  PRED_LVR = A_LVR/VLVR
  PRED_GUT = A_GUT/VGUT
  PRED_SPLR = A_SPLR/(0.25*VSPL)
  PRED_SPLS = A_SPLS/(0.75*VSPL)
  PRED_KID = A_KID/VKID
  PRED_HRT = A_HRT/VHRT
  PRED_MSC = A_MSC/VMSC
  PRED_BOD = A_BOD/VBOD
  
  IPRE_PA = PRED_PA*EPS
  IPRE_LNG = PRED_LNG*EPS
  IPRE_BRA = PRED_BRA*EPS
  IPRE_LVR = PRED_LVR*EPS
  IPRE_GUT = PRED_GUT*EPS
  IPRE_SPLR = PRED_SPLR*EPS
  IPRE_SPLS = PRED_SPLS*EPS
  IPRE_KID = PRED_KID*EPS
  IPRE_HRT = PRED_HRT*EPS
  IPRE_MSC = PRED_MSC*EPS
  IPRE_BOD = PRED_BOD*EPS
  
  Y = IPRE_PA

$THETA
 ; regional blood flow
  (0.4613) FIX  ; QBRA_STD 
  (2.251) FIX  ; QLVR_STD
  (1.5) FIX  ; QGUT_STD
  (0.1573) FIX  ; QSPL_STD
  (0.181) FIX  ; PSSPL_STD
  (1.272) FIX  ; QKID_STD
  (0.9227) FIX  ; QHRT_STD
  (2.223) FIX  ; QMSC_STD
  
 ; tissue mass balance (ml)
  (1.225) FIX  ; VMIX_STD
  (0.175) FIX  ; VLNG_STD
  (0.425) FIX  ; VBRA_STD
  (1.375) FIX  ; VLVR_STD
  (1.425) FIX  ; VGUT_STD
  (0.0875) FIX  ; VSPL_STD
  (0.425) FIX  ; VKID_STD
  (0.125) FIX  ; VHRT_STD
  (9.6) FIX  ; VMSC_STD
  
 ; renal physiology
  (0.33) FIX  ; GFR_STD
  (0.5) FIX  ; PSDIFF_STD
  (0.00111) FIX  ; KURINE_STD
  
 ; drug related parameters
  (0.001,0.6,0.999)  ; FU
  (6) FIX  ; VMAX
  (0.01) FIX  ; KM
  (0.001,0.03,99)  ; KA_STD
  (0.001,0.6,0.999)  ; FUTLNG
  (0.001,0.6,0.999)  ; FUTBRA
  (0.001,0.6,0.999)  ; FUTLVR
  (0.001,0.6,0.999)  ; FUTGUT
  (0.001,0.6,0.999)  ; FUTSPL
  (0.001,0.6,0.999)  ; FUTKID
  (0.001,0.6,0.999)  ; FUTHRT
  (0.001,0.6,0.999)  ; FUTMSC
  (0.001,0.6,0.999)  ; FUTBOD

$OMEGA 
  0 FIX ; PPVFUT

$SIGMA 
  1 ; EPS

$EST PRINT=1 METHOD=0 MAXEVALS=9999 POSTHOC NOABORT NSIG=3 SIGL=9
$COV

$TABLE ID TIME IPRE_PA IPRE_LNG IPRE_BRA IPRE_LVR IPRE_GUT IPRE_SPLR IPRE_SPLS 
IPRE_KID IPRE_HRT IPRE_MSC IPRE_BOD PRED_PA PRED_LNG PRED_BRA PRED_LVR PRED_GUT
PRED_SPLR PRED_SPLS PRED_KID PRED_HRT PRED_MSC PRED_BOD
WT TADNOM DOSEMGKG DOSE ROUTE
ONEHEADER NOPRINT FILE=lenalidomide.fit
