if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/.ihep.ac.cn/bes3/offline/ExternalLib/contrib/CMT/v1r20p20081118
endif
source ${CMTROOT}/mgr/setup.csh
set tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=KSKPiPiAlg -version=KSKPiPiAlg-00-00-01 -path=/besfs/users/caih/663p01/mywork $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}

