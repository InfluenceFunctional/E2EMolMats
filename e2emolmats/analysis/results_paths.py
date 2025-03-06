acridine_melt_paths = [
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface8/', # dev run
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface9/',  # convergence tests
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface10/',  # more convergence tests
    r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface11/',  # more convergence tests

    # old
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface1/',
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface2/', # failed
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface3/', # failed
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface5/',
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface6/', # something really weird happened here

    # old acridine ff
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface14/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface15/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface16_1/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface16_2/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface16_3/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface16_4/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface17_1/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface17_3/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface17_4/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface18/',
    #r'D:\crystal_datasets\acridine_melt_interface19/', # anthracene
    #r'D:\crystal_datasets\acridine_melt_interface20/'  # 2,7-DHN
]
acridine_scan_paths = [
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_interface_scan2/',
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_interface_scan3/', # first successful scan batch, with some refreezing
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_interface_scan4/', # single test
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_interface_scan5/',  # shorter test to compare new thermostat
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_interface_scan6/',  # different langevin dampings
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_interface_scan7/',  # different langevin dampings
    r'D:\crystal_datasets\acridine_w_new_ff/acridine_interface_scan8/',  # 2&4 melts
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_interface_scan9/',  # 3,6,7,8,9 melts folowing run 8
    r'D:\crystal_datasets\acridine_w_new_ff/acridine_interface_scan10/',  # 4 with different params

]
acridine_cluster_paths = [
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_cluster1/',  # dev run
    r'D:\crystal_datasets\acridine_w_new_ff/acridine_cluster2/',  # convergence run
    r'D:\crystal_datasets\acridine_w_new_ff/acridine_cluster3/',  # more convergence run

    # old acridine ff
    # r'D:\crystal_datasets\acridine_cluster4/',
    # r'D:\crystal_datasets\acridine_cluster5/',
    # r'D:\crystal_datasets\acridine_cluster6/',
    # r'D:\crystal_datasets\acridine_cluster7/',
    # r'D:\crystal_datasets\acridine_cluster8/',
    # r'D:\crystal_datasets\acridine_cluster9/',
    # r'D:\crystal_datasets\acridine_cluster10/',
    # r'D:\crystal_datasets\acridine_cluster11/',
    # r'D:\crystal_datasets\acridine_cluster12/',
    # r'D:\crystal_datasets\acridine_cluster13/',  # form 9 melt fix
    # r'D:\crystal_datasets\acridine_cluster14/',  # long runs
    # r'D:\crystal_datasets\acridine_cluster15/',  # init 27DHN runs
    # r'D:\crystal_datasets\acridine_cluster15_retest/',  # trying to rerun 15, where many runs failed

]
acridine_latent_paths = [
    # old acridine ff
    r'D:\crystal_datasets\acridine_w_old_ff/acridine_latents_battery1/',
    r'D:\crystal_datasets\acridine_w_old_ff/acridine_latents_battery2/',
]
acridine_cp_paths = [
    # old - Daisuke
    'D:\crystal_datasets\daisuke_cp_runs'
]
acridine_cp2_paths = [
    # r'D:\crystal_datasets\acridine_w_new_ff\acridine_cp3',  # dev run
    r'D:\crystal_datasets\acridine_w_new_ff\acridine_cp4',  # convergence test
    r'D:\crystal_datasets\acridine_w_new_ff\acridine_cp5',  # production runs

    # old runs
    #r'D:\crystal_datasets\acridine_w_new_ff\acridine_cp1',
    # r'D:\crystal_datasets\acridine_w_new_ff\acridine_cp2',

    ##old acridine ff
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_cp1',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_cp2',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_cp3',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_latents_battery1/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_latents_battery2/',
]
acridine_lattice_energy_paths = [
    r'D:\crystal_datasets\acridine_w_new_ff\acridine_lattice_energy1',  # gas phases
    r'D:\crystal_datasets\acridine_w_new_ff\acridine_lattice_energy2',  # solids
    r'D:\crystal_datasets\acridine_w_new_ff\acridine_lattice_energy3',  # gas phases
    r'D:\crystal_datasets\acridine_w_new_ff\acridine_lattice_energy4',  # gas phases

]
