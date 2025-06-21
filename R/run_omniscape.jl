#!/usr/bin/julia
using Omniscape
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--ini", "-i"
            help = "ini file path"
            required = true
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
ini_file = parsed_args["ini"]
run_omniscape(ini_file)
