%Copyright (C) 2022 by Frida Viset

run('tools\load_tools.m');
run('ModelShipExperiments\main.m');
run('ModelShipExperiments\VaryOdometryNoise.m');
run('ModelShipExperiments\RepeatedExperiments.m');
run('EKFSLAMOpenShoe\main.m');
run('Simulations\Vary_lengthscales.m');
run('Simulations\Vary_pred_pos_error.m');



