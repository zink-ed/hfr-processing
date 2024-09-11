from hfradarpy.radials import Radial

# data_path = '/home/cqiao/HFR_proc/hfradarpy/examples/data/'
# data_path = '/home/cqiao/media/disk/WFSC_MARA_MEAS_2024_04_Apr/'
# radial_file = data_path + 'radials/ruv/SEAB/' + 'RDLi_SEAB_2019_01_01_0200.ruv'
# data_path = '/home/cqiao/HFR_proc/'
raw_file = 'RDLm_MARA_2024_04_18_2200.ruv'
r0_file = 'RDLm_MARA_2024_04_18_2100.ruv'
radial_file = './radials_raw/' + raw_file
radial_file

# radial_file_qc = data_path + 'radials_qc/ruv/SEAB/' + 'RDLi_SEAB_2019_01_01_0200.ruv'
radial_file_qc = './radials_qc/' + raw_file
radial_file_qc

r = Radial(radial_file)

r.file_name

r.file_path

r.file_type()

r.metadata

r.data

r.diagnostics_hardware

r.diagnostics_radial

r.initialize_qc()

r.qc_qartod_syntax()
r.data.head()

r.qc_qartod_maximum_velocity(max_speed=250, high_speed=150) #cm/s
r.data.head()

r.qc_qartod_valid_location()
r.data.head()

r.qc_qartod_radial_count(low_count=140, min_count=50)
r.data.head()

r.qc_qartod_spatial_median(smed_range_cell_limit=2.1, smed_angular_limit=10, smed_current_difference=30)
r.data.head()

# r0 set so we can demonstrate temporal gradient test
# r0 = data_path + 'radials/ruv/SEAB/' + 'RDLi_SEAB_2019_01_01_0100.ruv'
r0 = './radials/' + r0_file

r.qc_qartod_temporal_gradient(r0, gradient_temp_fail=54, gradient_temp_warn=36)
r.data.head()

r.qc_qartod_avg_radial_bearing(reference_bearing=151, warning_threshold=15, failure_threshold=30)
r.data.head()

r.qc_qartod_primary_flag()
r.data.head()

r.to_ruv(radial_file_qc)
