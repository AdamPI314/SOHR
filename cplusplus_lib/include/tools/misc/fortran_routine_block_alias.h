#ifndef __FORTRAN_ROUTINE_BLOCK_ALIAS_H_
#define __FORTRAN_ROUTINE_BLOCK_ALIAS_H_

#include "global_macros.h"

#ifdef __WINDOWS_
#define ckstrt CKSTRT
#define chemkin CHEMKIN
#define ckstore CKSTORE
#define lsodestore LSODESTORE

#define chemkininitialize CHEMKININITIALIZE
#define ckrp CKRP
#define cknu CKNU
#define ckitr CKITR
#define ckncf CKNCF
#define ckrhoy CKRHOY
#define ckxty CKXTY
#define ckytx CKYTX
#define ckctx CKCTX
#define ckcty CKCTY
#define ckxtcp CKXTCP
#define ckytcr CKYTCR
#define cppdlsodev CPPDLSODEV
#define cppdlsodav CPPDLSODAV
#define cppdlsodevt CPPDLSODEVT
#define cppdlsodavt CPPDLSODAVT
#define cppdlsodep CPPDLSODEP
#define cppdlsodap CPPDLSODAP
#define cppdlsodept CPPDLSODEPT
#define cppdlsodapt CPPDLSODAPT
#define cppdlsodest CPPDLSODEST
#define cppdlsodast CPPDLSODAST
#define cppdlsodastcc1 CPPDLSODASTCC1
#define cppdlsodastcc2 CPPDLSODASTCC2
#define ckkfkr CKKFKR
#define ckkfkrsr CKKFKRSR
#define ckpy CKPY
#define ckcdyr CKCDYR
#define ckkfrt CKKFRT
#define ckraex CKRAEX
#define ckindx CKINDX
#define ckabe CKABE
#define ckcdc CKCDC
#define calculatetdotv CALCULATETDOTV

#endif



#ifdef __LINUX_

#define ckstrt ckstrt_
#define chemkin chemkin_
#define ckstore ckstore_
#define lsodestore lsodestore_

#define chemkininitialize chemkininitialize_
#define ckrp ckrp_
#define cknu cknu_
#define ckitr ckitr_
#define ckncf ckncf_
#define ckrhoy ckrhoy_
#define ckxty ckxty_
#define ckytx ckytx_
#define ckctx ckctx_
#define ckcty ckcty_
#define ckxtcp ckxtcp_
#define ckytcr ckytcr_
#define cppdlsodev cppdlsodev_
#define cppdlsodav cppdlsodav_
#define cppdlsodevt cppdlsodevt_
#define cppdlsodavt cppdlsodavt_
#define cppdlsodep cppdlsodep_
#define cppdlsodap cppdlsodap_
#define cppdlsodept cppdlsodept_
#define cppdlsodapt cppdlsodapt_
#define cppdlsodest cppdlsodest_
#define cppdlsodast cppdlsodast_
#define cppdlsodastcc1 cppdlsodastcc1_
#define cppdlsodastcc2 cppdlsodastcc2_
#define ckkfkr ckkfkr_
#define ckkfkrsr ckkfkrsr_
#define ckpy ckpy_
#define ckcdyr ckcdyr_
#define ckkfrt ckkfrt_
#define ckraex ckraex_
#define ckindx ckindx_
#define ckabe ckabe_
#define ckcdc ckcdc_
#define calculatetdotv calculatetdotv_

#endif


#endif
