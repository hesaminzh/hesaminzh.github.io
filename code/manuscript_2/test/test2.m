
I_panels = 2;
d = 2;
RIS_ref = [50, 0];

center_offset  = -(I_panels-1)/2;
panel_idx = (0:I_panels-1).';

build_panel_positions = @(dval) [RIS_ref(1) * ones(I_panels,1), RIS_ref(2) + (center_offset + panel_idx) * dval];

for i = 1:I_panels
    panel_positions(i,:) = [RIS_ref(1), RIS_ref(2) + (center_offset + (i-1)) * d];
end
