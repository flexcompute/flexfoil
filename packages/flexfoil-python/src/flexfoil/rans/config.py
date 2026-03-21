"""Flow360 case configuration builder for pseudo-2D airfoil RANS."""

from __future__ import annotations


def build_case_config(
    *,
    alpha: float,
    Re: float,
    mach: float,
    chord: float = 1.0,
    span: float = 0.01,
    temperature: float = 288.15,
    turbulence_model: str = "SpalartAllmaras",
    max_steps: int = 5000,
    cfl_initial: float = 5.0,
    cfl_final: float = 200.0,
    cfl_ramp_steps: int = 2000,
    order_of_accuracy: int = 2,
) -> dict:
    """Build a Flow360 case JSON configuration for pseudo-2D RANS.

    Parameters
    ----------
    alpha : float
        Angle of attack in degrees.
    Re : float
        Reynolds number based on chord.
    mach : float
        Freestream Mach number.
    chord : float
        Chord length (default 1.0, nondimensional).
    span : float
        Spanwise extent of the pseudo-3D mesh (default 0.01).
    temperature : float
        Freestream temperature in Kelvin (default 288.15 K = ISA sea level).
    turbulence_model : str
        'SpalartAllmaras' or 'kOmegaSST'.
    max_steps : int
        Maximum pseudo-time steps for steady convergence.
    cfl_initial, cfl_final, cfl_ramp_steps : float
        CFL number ramping schedule.
    order_of_accuracy : int
        Spatial order of accuracy (1 or 2).

    Returns
    -------
    dict
        Flow360 case JSON configuration.
    """
    ref_area = chord * span  # wetted area for 2D coefficient normalization

    config = {
        "geometry": {
            "refArea": ref_area,
            "momentCenter": [chord * 0.25, 0.0, 0.0],
            "momentLength": [chord, chord, chord],
        },
        "freestream": {
            "Mach": mach,
            "Reynolds": Re,
            "alphaAngle": alpha,
            "betaAngle": 0.0,
            "Temperature": temperature,
        },
        # UGRID boundary tags are integers: 1=wall, 2=farfield, 3=sym_z0, 4=sym_z1
        "boundaries": {
            "1": {"type": "NoSlipWall"},
            "2": {"type": "Freestream"},
            "3": {"type": "SlipWall"},
            "4": {"type": "SlipWall"},
        },
        "navierStokesSolver": {
            "absoluteTolerance": 1e-10,
            "linearIterations": 35,
            "kappaMUSCL": -1.0,
            "orderOfAccuracy": order_of_accuracy,
        },
        "turbulenceModelSolver": {
            "modelType": turbulence_model,
            "absoluteTolerance": 1e-8,
            "linearIterations": 25,
            "orderOfAccuracy": order_of_accuracy,
        },
        "timeStepping": {
            "maxPhysicalSteps": 1,
            "maxPseudoSteps": max_steps,
            "timeStepSize": "inf",
            "CFL": {
                "initial": cfl_initial,
                "final": cfl_final,
                "rampSteps": cfl_ramp_steps,
            },
        },
        "surfaceOutput": {
            "outputFormat": "paraview",
            "animationFrequency": -1,
            "surfaces": {
                "1": {
                    "outputFields": ["Cp", "Cf", "CfVec", "yPlus"],
                },
            },
        },
        "volumeOutput": {
            "outputFormat": "paraview",
            "animationFrequency": -1,
            "outputFields": ["primitiveVars", "Mach", "Cp"],
        },
    }

    return config
