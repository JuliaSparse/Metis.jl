using BinDeps

@BinDeps.setup

libmetis = library_dependency("libmetis", aliases=["libmetis5"])

#provides(AptGet,"libmetis5",libmetis)
#provides(Yum,"metis-5.1.0",libmetis)

provides(Sources,
         URI("http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz"),
         libmetis)

metisbuilddir = BinDeps.builddir(libmetis)
println(libmetis)
println(libdir(libmetis))
println(metisbuilddir)

provides(BuildProcess,
	(@build_steps begin
            println("In buildsteps")
	    GetSources(libmetis)
	    println("past GetSources")
	    CreateDirectory(metisbuilddir)
            println("past CreateDirectory")
	    @build_steps begin
		ChangeDirectory(metisbuilddir)
	    end
	end),libmetis)

@BinDeps.install
