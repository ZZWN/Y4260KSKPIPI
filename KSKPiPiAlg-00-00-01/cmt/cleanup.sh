if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/.ihep.ac.cn/bes3/offline/ExternalLib/contrib/CMT/v1r20p20081118; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then tempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=KSKPiPiAlg -version=KSKPiPiAlg-00-00-01 -path=/besfs/users/caih/663p01/mywork $* >${tempfile}; . ${tempfile}
/bin/rm -f ${tempfile}

