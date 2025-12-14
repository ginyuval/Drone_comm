function a = steering_vec_ula(thetaDeg, gl_params)
%STEERING_VEC_ULA ULA steering vector using sin(theta) convention.
    thetaRad = deg2rad(thetaDeg);
    m = (0:gl_params.numElements-1).';
    a = exp(1j * 2*pi * gl_params.d / gl_params.lambda * m * sin(thetaRad));
end
