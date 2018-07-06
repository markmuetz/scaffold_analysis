from omnium.converter import FF2NC_Converter


class ExptConverter(FF2NC_Converter):
    analysis_name = 'expt_converter'
    single_file = True
    input_dir = 'work/20000101T0000Z/{expt}_atmos'
    input_filename = '{input_dir}/atmos.pp1'
    output_dir = 'work/20000101T0000Z/{expt}_atmos'
    output_filenames = ['{output_dir}/atmos.pp1.nc']

    force = True
    delete = True
