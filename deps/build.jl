using BinDeps

@BinDeps.setup

libmetis = library_dependency("libmetis", aliases=["libmetis5"])

@windows_only begin
  using WinRPM
  provides(WinRPM.RPM, "metis", libmetis, os = :Windows)
end

@osx_only begin
    if Pkg.installed("Homebrew") === nothing
        error("Homebrew package not installed, please run Pkg.add(\"Homebrew\")")  end
    using Homebrew
  provides(Homebrew.HB, "metis", libmetis, os = :Darwin)
end

provides(AptGet, "libmetis5", libmetis)

provides(Yum, "metis-5.1.0", libmetis)

provides(Sources, URI("http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz"), libmetis)

metisname = "metis-5.1.0"

libdir = BinDeps.libdir(libmetis)
downloadsdir = BinDeps.downloadsdir(libmetis)
builddir = BinDeps.builddir(libmetis)

provides(SimpleBuild,
         (@build_steps begin
             GetSources(libmetis)
             CreateDirectory(builddir)
             (@build_steps begin
	         ChangeDirectory(builddir)
	         (@build_steps begin
	             `make config shared=1 prefix=$libdir/..`
	             `make`
	             `make install`
                 end)
	     end)
         end), [libmetis], os=:Unix)

@BinDeps.install [:libmetis => :libmetis]
