from rsqsim_api.fault.multifault import RsqSimMultiFault, RsqSimSegment

mesh = "slip_distribution.vtk"
seg = RsqSimSegment.from_vtk(mesh)
seg.name = "variable_dip"
seg.to_rsqsim_fault_file("variable_dip.flt", mm_yr=True)

mesh = "slip_distribution_ndc.vtk"
seg = RsqSimSegment.from_vtk(mesh)
seg.name = "no_dip_change"
seg.to_rsqsim_fault_file("no_dip_change.flt", mm_yr=True)