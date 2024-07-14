## Queueing 
Checking the queue with
```
squeue
```
example output:
```
JOBID PARTITION     NAME     USER     ST   TIME  NODES NODELIST(REASON)
4     h3            jet_4_mc omry_mag PD   0:00      1 (BeginTime)
```

- **JOBID**: `4` - This is the unique identifier for your job.
- **PARTITION**: `h3` - The partition you requested for the job.
- **NAME**: `jet_4_mc` - The name of your job.
- **USER**: `omry_mag` - Your username.
- **ST**: `PD` - The status of the job. `PD` stands for "Pending".
- **TIME**: `0:00` - The amount of time the job has been running. Since it is still pending, this is `0:00`.
- **NODES**: `1` - The number of nodes allocated for the job.
- **NODELIST(REASON)**: `(BeginTime)` - This indicates why the job is pending, here, waiting for begin time.
### Status (ST) options
- **PD (PENDING)**: The job is waiting for resource allocation. This is typically due to job dependencies, resource availability, or scheduling policies.
- **R (RUNNING)**: The job is currently being executed.
- **CG (COMPLETING)**: The job is in the process of completing. The job script has finished, but job cleanup is still in progress.
- **CD (COMPLETED)**: The job has finished successfully.
- **F (FAILED)**: The job terminated with a non-zero exit code or was otherwise unsuccessful.
- **TO (TIMEOUT)**: The job exceeded its time limit and was terminated by the system.
- **CA (CANCELLED)**: The job was cancelled by the user or system administrator.
- **NF (NODE_FAIL)**: The job terminated due to a failure of one or more allocated nodes.
- **SE (SPECIAL_EXIT)**: The job terminated with a special exit state, often used for jobs that are part of a larger job array.
- **RV (REVOKED)**: The job was revoked, typically due to preemption or a system reservation.
- **OOM (OUT_OF_MEMORY)**: The job was terminated because it exceeded the memory limit.
### When is my job scheduled to begin?

write in the following command with the given job ID. In this case, it is 4:

```
scontrol show job 4
```

output:
```
   UserId=omry_magen_gmail_com(1536403189) GroupId=omry_magen_gmail_com(1536403189) MCS_label=N/A
   Priority=1 Nice=0 Account=(null) QOS=normal
   JobState=PENDING Reason=BeginTime Dependency=(null)
   Requeue=1 Restarts=16 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:00 TimeLimit=00:01:00 TimeMin=N/A
   SubmitTime=2024-07-14T12:46:41 EligibleTime=2024-07-14T12:48:42
   AccrueTime=2024-07-14T12:48:42
   StartTime=2024-07-14T12:48:42 EndTime=2024-07-14T12:49:42 Deadline=N/A
   SuspendTime=None SecsPreSuspend=0 LastSchedEval=2024-07-14T12:46:37 Scheduler=Main
   Partition=h3 AllocNode:Sid=hpcomry-slurm-login-001:2343
   ReqNodeList=(null) ExcNodeList=(null)
   NodeList=
   NumNodes=1-1 NumCPUs=16 NumTasks=16 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   ReqTRES=cpu=16,mem=63536M,node=1,billing=16
   AllocTRES=cpu=88,mem=349448M,node=1,billing=88
   Socks/Node=* NtasksPerN:B:S:C=16:0:*:* CoreSpec=*
   MinCPUsNode=16 MinMemoryCPU=3971M MinTmpDiskNode=0
   Features=(null) DelayBoot=00:00:00
   OverSubscribe=NO Contiguous=0 Licenses=(null) Network=(null)
   Command=/home/omry_magen_gmail_com/workdir/16_tasks/jet.slurm
   WorkDir=/home/omry_magen_gmail_com/workdir/16_tasks
   AdminComment=GCP Error: No eligible zone could be found in this region for given properties 
   StdErr=/home/omry_magen_gmail_com/workdir/16_tasks/slurm-4.out
   StdIn=/dev/null
   StdOut=/home/omry_magen_gmail_com/workdir/16_tasks/slurm-4.out
   Power=
   TresPerTask=cpu:1
```

This summarises the run script that we launched before. Most importantly, we can see that the job was scheduled to start 2 minutes after it was submitted.


```
scontrol show node hpcomry-h3nodeset-2
```