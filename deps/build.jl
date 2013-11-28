using BinDeps

@BinDeps.setup

deps = [
        metis = library_dependency("metis", aliases=["libmetis","libmetis5"])
        ]

#provides(AptGet,"libmetis5",metis)


# version of metis package to use
const metisver = "5.1.0"
const metisdir = string("metis-", metisver)
const metisdownload = string(metisdir, ".tar.gz")
const metisurl = string("http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/", metisdownload)
provides(Sources,URI(metisurl),metis)


@BinDeps.install [:metis => :_jl_metis]

