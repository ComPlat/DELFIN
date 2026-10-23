"""The OpenMPI environment DELFIN runs parallel ORCA under, in one place.

Inside a SLURM allocation OpenMPI takes its slot count from the resource
manager: ``--ntasks=1`` is one slot, however many CPUs the task has, and a
parallel ORCA asking for eight processes on sixteen idle cores dies with
"There are not enough slots available in the system" -- measured on this
cluster (2026-09-26), it killed every band and every OptTS the dashboard's
saddle search started, in seconds, on a machine with nothing else to do.

DELFIN's job runners have carried these settings for exactly that reason
(:mod:`delfin.dashboard.local_runner` and the submit templates beside it).
They live here now so that the one other place that starts a parallel ORCA
-- the editor's saddle search -- sets the same environment from the same
list rather than a copy of it, because two copies of a list of OpenMPI
knobs is one copy more than anyone keeps current.

The values are what the job runners set, unchanged: the ob1 point-to-point
layer with the self, tcp and vader transports, no hcoll, no hwloc binding
(saddle searches share a node with the dashboard itself), core mapping with
oversubscribe on -- the setting that answers the one-task allocation -- and
a yield when idle, because an oversubscribed run that spins helps nobody.
"""

#: Set on the environment of every parallel ORCA DELFIN starts, unless the
#: user set the same variable themselves: ``setdefault`` semantics, so a
#: site's own MCA tuning wins over ours.
OMPI_ENVIRONMENT = {
    'OMPI_MCA_pml': 'ob1',
    'OMPI_MCA_btl': 'self,tcp,vader',
    'OMPI_MCA_mpi_show_mca_params_file': '0',
    'OMPI_MCA_mpi_yield_when_idle': '1',
    'OMPI_MCA_coll_hcoll_enable': '0',
    'OMPI_MCA_hwloc_base_binding_policy': 'none',
    'OMPI_MCA_rmaps_base_mapping_policy': 'core',
    'OMPI_MCA_rmaps_base_oversubscribe': 'true',
}

#: Taken off again, because a cluster-wide default for either of these has
#: broken more runs than it fixed: the include list names interfaces this
#: cluster does not have, and the exclude list names one the transport does.
OMPI_ENVIRONMENT_REMOVED = (
    'OMPI_MCA_btl_tcp_if_include',
    'OMPI_MCA_btl_tcp_if_exclude',
)
