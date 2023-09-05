using Sunny


if Sys.iswindows()
    SunnyAssist_path = "C:\\Users\\vdp\\Dropbox (ORNL)\\Sunny\\SunnyAssist\\"
else
    SunnyAssist_path = "/home/vdp/Dropbox (ORNL)/Sunny/SunnyAssist/"
end

include(joinpath(SunnyAssist_path, "SunnyAssist.jl"))

if Sys.iswindows()
    path = "C:\\Users\\vdp\\Dropbox (ORNL)\\Sunny\\SunnyAssist\\TBFO_example\\"
else
    path = "/home/vdp/Dropbox (ORNL)/Sunny/SunnyAssist/TBFO_example/"
end

RegionA = load_nxs(joinpath(path, "RegionA.nxs"))
RegionB = load_nxs(joinpath(path, "RegionB.nxs"))
RegionC = load_nxs(joinpath(path, "RegionC.nxs"))

energy_scan_H0p5_K1p5_Lm0p2 = load_nxs(joinpath(path, "energy_scan_H0p5_K1p5_Lm0p2.nxs"))
energy_scan_H0p5_K1p5_Lm0p2 = load_nxs(joinpath(path, "energy_scan_H0p5_K1p5_Lm1p6.nxs"))
energy_scan_H0p5_K1p5_Lm0p2 = load_nxs(joinpath(path, "energy_scan_H0p5_K1p5_Lm2.nxs"))
energy_scan_H0p5_K1p5_Lm0p2 = load_nxs(joinpath(path, "energy_scan_H1p3_K0p5_Lm2.nxs"))
energy_scan_H0p5_K1p5_Lm0p2 = load_nxs(joinpath(path, "energy_scan_H1p5_K0p5_Lm2.nxs"))
energy_scan_H0p5_K1p5_Lm0p2 = load_nxs(joinpath(path, "energy_scan_H1p9_K0p5_Lm1.nxs"))
energy_scan_H0p5_K1p5_Lm0p2 = load_nxs(joinpath(path, "energy_scan_Hm0p5_K0p5_Lm1.nxs"))
energy_scan_H0p5_K1p5_Lm0p2 = load_nxs(joinpath(path, "energy_scan_Hm0p5_K0p7_Lm1.nxs"))


print(RegionA[1])