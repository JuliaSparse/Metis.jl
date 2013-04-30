using BinDeps

# version of metis package to use
const metisver = "5.1.0"
const metisdir = string("metis-", metisver)
const metisdownload = string(metisdir, ".tar.gz")
const metisurl = string("http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/", metisdownload)

@unix_only begin
    prefix = joinpath(Pkg.dir(),"Metis","deps","usr")
    if !isfile(metisdownload)
        global cmkverstring=""
        global cmkver=Array(Int,4)
        try
            cmkverstring = split(readall(`cmake --version`))[3]
            cmkver = map(int, split(cmkverstring, '.'))
        catch
            error("Please install cmake version 2.8 or later")
        end
        if cmkver[1] < 2 || cmkver[2] < 8
            error("cmake is version $cmkverstring, version 2.8 or later is required")
        end
        run(download_cmd(metisurl, metisdownload))
        run(unpack_cmd(metisdownload, "."))
        cd(metisdir) do
            run(`make config shared=1 prefix=$prefix`)
            run(`make install`)
            run(`make clean`)
        end
    end
end # unix_only
