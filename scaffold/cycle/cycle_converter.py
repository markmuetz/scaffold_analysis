from omnium.converter import FF2NC_Converter


class CycleConverter(FF2NC_Converter):
    analysis_name = 'cycle_converter'
    single_file = True
    input_dir = 'share/data/history/{expt}'
    input_filename_glob = '{input_dir}/atmos.???.pp1'
    output_dir = 'work/20000101T0000Z/{expt}_atmos'
    output_filenames = ['{output_dir}/atmos.pp1.nc']

    uses_runid = True
    runid_pattern = 'atmos.(?P<runid>\d{3}).pp1'

    force = True
    delete = False
