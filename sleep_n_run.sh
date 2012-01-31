while true;do
  while qstat | grep 'phillict' | grep 'parallel_sa ' | egrep ' [RQ] ' | grep janus-normal;do 
    echo "Waiting for job to complete"
    jid=$(qstat | grep 'phillict' | grep 'parallel_sa ' | egrep ' [RQ] ' | grep janus-normal | cut -f 1 -d ' ')
    echo "================================= CHECKJOB ==================================="
    checkjob $jid
    echo "================================= SHOWSTART =================================="
    showstart $jid
    id=$(echo $jid | cut -f 1 -d '.')
    if [ -f parallel_sa.o$id ];then
      echo "==================================== LOG ====================================="
      tail -n 10 "parallel_sa.o$id"
      echo "=============================================================================="
    fi
    sleep 5m
  done
  echo "Running new job"
  qsub run.sh
  sleep 1h
done
