function feature_hslice = return_h_slice(param, ii_h, ii_knot, ii_dofs)

    if isempty(ii_h)
        ii_h = 1:size(param.q_h_array, 3);
    end

    feature_hslice.q =        param.q_h_array   (:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.dq =       param.dq_h_array  (:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.ddq =      param.ddq_h_array (:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.dddq =     param.dddq_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.x =        param.x_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.dx =       param.dx_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.ddx =      param.ddx_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.dddx =     param.dddx_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.tau =      param.tau_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.dtau =     param.dtau_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.ddtau =    param.ddtau_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.cop =      param.cop_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.dcop =     param.dcop_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.ddcop =    param.ddcop_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.com =      param.com_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.dcom =     param.dcom_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.ddcom =    param.ddcom_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.ep =       param.ep_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.ek =       param.ek_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.geo =      param.geo_h_array(:, :, ii_h, ii_knot, ii_dofs);
    feature_hslice.en =       param.en_h_array(:, :, ii_h, ii_knot, ii_dofs);
