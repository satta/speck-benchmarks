#!/bin/bash

# from http://stackoverflow.com/questions/1080461/
#      /peak-memory-measurement-of-long-running-process-in-linux

ppid=$$

# VmPeak: Peak virtual memory size.
# VmSize: Virtual memory size.
# VmLck: Locked memory size.
# VmHWM: Peak resident set size ("high water mark").
# VmRSS: Resident set size.
# VmData, VmStk, VmExe: Size of data, stack, and text segments.
# VmLib: Shared library code size.
# VmPTE: Page table entries size (since Linux 2.6.10).

maxVmPeak=0
maxVmSize=0
maxVmLck=0
maxVmHWM=0
maxVmRSS=0
maxVmData=0
maxVmStk=0
maxVmExe=0
maxVmLib=0
maxVmPTE=0

/usr/bin/time --format="time:%e" $@ 2>&1 &
tpid=`pgrep -P ${ppid} -n -f time`
if [[ ${tpid} -ne "" ]]; then
  pid=`pgrep -P ${tpid} -n -f $1` # $! may work here but not later
fi

while [[ ${tpid} -ne "" ]]; do
  if [[ ${pid} -ne "" ]]; then
    VmPeak=`cat /proc/${pid}/status | grep VmPeak | awk '{print $2}'`
    if [[ ${VmPeak} -gt ${maxVmPeak} ]]; then
      maxVmPeak=${VmPeak}
    fi
    VmSize=`cat /proc/${pid}/status | grep VmSize | awk '{print $2}'`
    if [[ ${VmSize} -gt ${maxVmSize} ]]; then
      maxVmSize=${VmSize}
    fi
    VmLck=`cat /proc/${pid}/status | grep VmLck | awk '{print $2}'`
    if [[ ${VmLck} -gt ${maxVmLck} ]]; then
      maxVmLck=${VmLck}
    fi
    VmHWM=`cat /proc/${pid}/status | grep VmHWM | awk '{print $2}'`
    if [[ ${VmHWM} -gt ${maxVmHWM} ]]; then
      maxVmHWM=${VmHWM}
    fi
    VmRSS=`cat /proc/${pid}/status | grep VmRSS | awk '{print $2}'`
    if [[ ${VmRSS} -gt ${maxVmRSS} ]]; then
      maxVmRSS=${VmRSS}
    fi
    VmData=`cat /proc/${pid}/status | grep VmData | awk '{print $2}'`
    if [[ ${VmData} -gt ${maxVmData} ]]; then
      maxVmData=${VmData}
    fi
    VmStk=`cat /proc/${pid}/status | grep VmStk | awk '{print $2}'`
    if [[ ${VmStk} -gt ${maxVmStk} ]]; then
      maxVmStk=${VmStk}
    fi
    VmExe=`cat /proc/${pid}/status | grep VmExe | awk '{print $2}'`
    if [[ ${VmExe} -gt ${maxVmExe} ]]; then
      maxVmExe=${VmExe}
    fi
    VmLib=`cat /proc/${pid}/status | grep VmLib | awk '{print $2}'`
    if [[ ${VmLib} -gt ${maxVmLib} ]]; then
      maxVmLib=${VmLib}
    fi
    VmPTE=`cat /proc/${pid}/status | grep VmPTE | awk '{print $2}'`
    if [[ ${VmPTE} -gt ${maxVmPTE} ]]; then
      maxVmPTE=${VmPTE}
    fi
  fi
  sleep 1
  savedtpid=${tpid}
  tpid=`pgrep -P ${ppid} -n -f time`
done
wait ${savedtpid} # don't wait, job is finished
exitstatus=$?   # catch the exit status of wait, the same of $@
echo "Memory usage for $@:"
echo "  VmPeak: ${maxVmPeak} kB"
echo "  VmSize: ${maxVmSize} kB"
echo "  VmLck:  ${maxVmLck} kB"
echo "  VmHWM:  ${maxVmHWM} kB"
echo "  VmRSS:  ${maxVmRSS} kB"
echo "  VmData: ${maxVmData} kB"
echo "  VmStk:  ${maxVmStk} kB"
echo "  VmExe:  ${maxVmExe} kB"
echo "  VmLib:  ${maxVmLib} kB"
echo "  VmPTE:  ${maxVmPTE} kB"
echo "Exit status: ${exitstatus}"

