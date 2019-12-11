function Q_ECI_RIC = eci2ric_dcm(r_eci,v_eci)
%ECI2RIC_DCM Earth Centered Inertial Frame to Radial, In-track, cross-track
%frame DCM 

% Unit Vectors, equal to rows of DCM
uhat = (r_eci./norm(r_eci))';
cr_prod = (cross(r_eci, v_eci))';
what = cr_prod./norm(cr_prod);
cr_prod = cross(what, uhat);
vhat = cr_prod./norm(cr_prod);

% Assemble dcm
Q_ECI_RIC = [uhat; what; vhat];
end

