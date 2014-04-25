using BinDeps

@BinDeps.setup

metis = library_dependency("libmetis", aliases=["libmetis5"])

@windows_only begin
  using WinRPM
  provides(WinRPM.RPM, "metis", metis, os = :Windows )
end

@osx_only begin
  using Homebrew
  provides( Homebrew.HB, "metis", metis, os = :Darwin )
end

provides( AptGet, "libmetis5", metis )
provides( Yum, "metis-5.1.0", metis )

julia_usrdir = normpath(JULIA_HOME*"/../") # This is a stopgap
libdirs = String["$(julia_usrdir)/lib"]
includedirs = String["$(julia_usrdir)/include"]

provides( Sources, URI("http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz"), metis )
## provides( BuildProcess,
##           Autotools(lib_dirs = libdirs,
##                     include_dirs = includedirs,
##                     configure_options = ["--shared=1"]),
##           metis )

@BinDeps.install [:metis => :metis]
