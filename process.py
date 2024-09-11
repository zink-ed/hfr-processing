from hfradarpy.radials import Radial, qc_radial_file
import glob
import os

folder = 'JEFF/'

# Path to radial directory
radial_dir = '/home/cqiao/HFR_proc/radials_raw/' + folder
save_dir = '/home/cqiao/HFR_proc/radials_qc/' + folder
clean_dir = '/home/cqiao/HFR_proc/radials_clean/' + folder

# Use glob to find radial files (*
files = sorted(glob.glob(os.path.join(radial_dir, '*.ruv')))

def run_tests(r, prev):

    # run high frequency radar qartod tests on open radial file

    qc_values = dict(
        qc_qartod_avg_radial_bearing=dict(reference_bearing=151, warning_threshold=15, failure_threshold=30),
        qc_qartod_radial_count=dict(min_count=75.0, low_count=225.0),
        qc_qartod_maximum_velocity=dict(max_speed=300.0, high_speed=100.0),
        qc_qartod_spatial_median=dict(smed_range_cell_limit=2.1, smed_angular_limit=10, smed_current_difference=30),
        qc_qartod_temporal_gradient=dict(gradient_temp_fail=32, gradient_temp_warn=25),
        qc_qartod_primary_flag=dict(include=['qc_qartod_syntax', 'qc_qartod_valid_location', 'qc_qartod_radial_count',
                                             'qc_qartod_maximum_velocity', 'qc_qartod_spatial_median'])
    )
    
    qc_radial_file(radial_file=r, qc_values=qc_values, export="radial", save_path=save_dir, clean=True, clean_path=clean_dir)
    
    #file_name = r.file_name[:-4:] + '_proc.ruv'

    # r.to_ruv(save_dir + file_name)
    
prev = files[0]

for f in files[1:]:
    r = Radial(f)
    print(r.file_name)
    run_tests(r, prev)
    prev = radial_dir + r.file_name
    print(prev)
