from omnium.converter import FF2NC_Converter


class ExptConverter(FF2NC_Converter):
    analysis_name = 'expt_converter'
    single_file = True
    input_dir = 'work/20000101T0000Z/{expt}_atmos'
    # A glob so that if it's not there it doesn't raise an error.
    input_filename_glob = '{input_dir}/atmos.pp3'
    output_dir = 'work/20000101T0000Z/{expt}_atmos'
    output_filenames = ['{output_dir}/atmos.pp3.nc']

    force = False
    delete = True
