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
#define ckcppdlsodev CKCPPDLSODEV
#define ckcppdlsodav CKCPPDLSODAV
#define ckcppdlsodevt CKCPPDLSODEVT
#define ckcppdlsodavt CKCPPDLSODAVT
#define ckcppdlsodep CKCPPDLSODEP
#define ckcppdlsodap CKCPPDLSODAP
#define ckcppdlsodept CKCPPDLSODEPT
#define ckcppdlsodapt CKCPPDLSODAPT
#define ckcppdlsodest CKCPPDLSODEST
#define ckcppdlsodast CKCPPDLSODAST
#define ckcppdlsodastcc1 CKCPPDLSODASTCC1
#define ckcppdlsodastcc2 CKCPPDLSODASTCC2
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
#define ckcppdlsodev ckcppdlsodev_
#define ckcppdlsodav ckcppdlsodav_
#define ckcppdlsodevt ckcppdlsodevt_
#define ckcppdlsodavt ckcppdlsodavt_
#define ckcppdlsodep ckcppdlsodep_
#define ckcppdlsodap ckcppdlsodap_
#define ckcppdlsodept ckcppdlsodept_
#define ckcppdlsodapt ckcppdlsodapt_
#define ckcppdlsodest ckcppdlsodest_
#define ckcppdlsodast ckcppdlsodast_
#define ckcppdlsodastcc1 ckcppdlsodastcc1_
#define ckcppdlsodastcc2 ckcppdlsodastcc2_
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


#ifdef __USE_CANTERA_

#define canterainitialize canterainitialize_
#define ctindx ctindx_
#define ctwyp ctwyp_
#define ctwyr ctwyr_
#define ctcdyr ctcdyr_
#define ctkfkr ctkfkr_

#define canteracppdlsodev canteracppdlsodev_
#define canteracppdlsodav canteracppdlsodav_


#endif

#endif
