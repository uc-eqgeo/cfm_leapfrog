from rsqsim_api.fault.multifault import RsqSimMultiFault, RsqSimSegment

mesh = "slip_distribution_truncated_north.vtk"
seg1 = RsqSimSegment.from_vtk(mesh, segment_number=0, fault_name="truncated_north")

mesh = "slip_distribution_truncated_south.vtk"
seg2 = RsqSimSegment.from_vtk(mesh, segment_number=1, fault_name="truncated_south")

combined = RsqSimMultiFault([seg1, seg2])
combined.write_rsqsim_input_file("truncated_faults.flt", mm_yr=True)