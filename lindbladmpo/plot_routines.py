# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""
Routines for basic plotting of simulation results, with both general and more specialized functions.
"""

from typing import Optional, Tuple, List, Union, Any, Sequence

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import numpy as np
from lindbladmpo.LindbladMPOSolver import LindbladMPOSolver


LINDBLADMPO_TEX_LABELS = {
    "tr_rho": "{\\rm tr}\\rho",
    "s_2": "S_2",
    "osee_center": "OSEE_{\\rm center}",
    "max_bond_dim": "{\\rm max bond-dim}",
    "duration_ms": "{\\rm duration(ms)}",
}
"""Labels in LaTex format for global output data, indexed by their data file entries.
	Note that key entries here are lower-case."""


LINDBLADMPO_LINE_STYLES = ["-", "--", ":", "-."]
"""Line styles for curves plotted using matplotlib"""


def prepare_time_data(
    parameters: dict,
    n_t_ticks=10,
    t_ticks_round=3,
    t_init: Optional[float] = None,
    t_final: Optional[float] = None,
) -> (np.ndarray, np.ndarray, np.ndarray, int):
    """
    Prepare the time variables used in plotting the results from one simulation.

    Args:
            parameters: A dictionary from which the basic time parameters are taken.
            n_t_ticks: The number of labeled major tick marks to generate for the time axis.
            t_ticks_round: The number of digits to round time axis tick labels to.
            t_init: An optional initial time - if None, it is taken from the parameters dict.
            t_final: An optional final time - if None, it is taken from the parameters dict.

    Returns:
            A tuple with the following four entries:
                    t_eval: An array of the times for which the solver parameters indicate output is to
                            be evaluated for (based on the `t_init`, `t_final`, and `tau` parameters.
                    t_tick_indices: The indices of the time axis tick marks.
                    t_tick_labels: The formatted labels of the time axis tick marks.
                    n_t_steps: The number of time steps
    """
    tau = parameters["tau"]
    if t_final is None:
        t_final = parameters["t_final"]
    if t_init is None:
        t_init = parameters.get("t_init", 0.0)
    # output_step = parameters.get('output_step', 1)
    # TODO handle output_step
    n_t_steps = int((t_final - t_init) / (tau * 1)) + 1
    t_tick_indices = np.arange(0, n_t_steps, int(n_t_steps / n_t_ticks))
    # TODO: Fix division by zero when n_t_ticks < n_t_steps above
    t_eval = np.arange(t_init, t_final, tau)
    if t_final not in t_eval:
        t_eval = np.concatenate((t_eval, [t_final]))
    if t_ticks_round == 0:
        t_tick_labels = np.asarray(t_init + t_tick_indices * tau, dtype=int)
    else:
        t_tick_labels = np.round(t_init + t_tick_indices * tau, t_ticks_round)
    return t_eval, t_tick_indices, t_tick_labels, n_t_steps


def prepare_curve_data(
    result: dict,
    s_output_type: str,
    s_obs_name: str,
    q_indices: Union[Tuple, Tuple[int]],
) -> ((list, list), str):
    """
    Prepare the data used for plotting one curve of simulation observables.

    Args:
            result: A dictionary from which the observables are taken.
            s_output_type: The type of output, used a key into the result dict, and also in formatting
                    the descriptive tex label of the data.
            s_obs_name: The name of the specific observable, used a key into the relevant observables
                    dict, and also in formatting the descriptive tex label of the data.
            q_indices: A tuple with the indices of the qubits identifying the observable to plot, or
                    an empty tuple if the observable is a global one.

    Returns:
            A tuple with the following two entries:
                    obs_data: A tuple of two lists, the first being the time points, the second being
                    the data.
                    s_tex_label: A formatted tex label for the data.
    """
    obs_dict = result[s_output_type]
    obs_data = None
    s_tex_label = ""
    s_obs_name = s_obs_name.lower()
    if obs_dict is not None:
        obs_data = obs_dict.get((s_obs_name, q_indices))
        if s_output_type == "obs-1q":
            s_tex_label = f"\\sigma^{s_obs_name}_{{{q_indices[0]}}}"
        elif s_output_type == "obs-2q":
            s_tex_label = f"\\sigma^{s_obs_name[0]}_{{{q_indices[0]}}} " \
                          f"\\sigma^{s_obs_name[1]}_{{{q_indices[1]}}}"
        elif s_output_type == "global":
            s_tex_label = LINDBLADMPO_TEX_LABELS[s_obs_name]
    return obs_data, s_tex_label


def prepare_2q_correlation_data(
    result: dict, s_obs_name: str, q_indices: Tuple[int]
) -> ((list, list), str):
    """
    Prepare the data used for plotting one curve of a two-qubit connected correlation function.
    The connected correlation is defined as a 2Q observable from which the product of the two
    corresponding 1Q observables is subtracted.

    Args:
            result: A dictionary from which the observables are taken.
            s_obs_name: The name of the specific observable, used a key into the relevant observables
                    dict, and also in formatting the descriptive tex label of the data.
            q_indices: A tuple with the indices of the qubits identifying the correlation function.

    Returns:
            A tuple with the following two entries:
                    obs_data: A tuple of two lists, the first being the time points, the second being
                            the data.
                    s_tex_label: A formatted tex label for the data - note that this string does not
                            indicate that the data is of 2Q connected correlation (the subtraction of
                            the 1Q product).
    """
    obs_1q_dict = result["obs-1q"]
    obs_2q_dict = result["obs-2q"]
    obs_data = None
    s_tex_label = ""
    s_obs_name = s_obs_name.lower()
    # Below we verify that all required 1Q and 2Q observables have complete data available.
    if obs_1q_dict is not None and obs_2q_dict is not None:
        obs_0 = obs_1q_dict.get((s_obs_name[0], (q_indices[0],)))
        obs_1 = obs_1q_dict.get((s_obs_name[1], (q_indices[1],)))
        obs_2 = obs_2q_dict.get((s_obs_name, q_indices))
        if (
            obs_0 is not None
            and obs_1 is not None
            and obs_2 is not None
            and len(obs_0[0]) == len(obs_1[0])
            and len(obs_0[0]) == len(obs_2[0])
        ):
            # The data all comes from one simulation, so we can safely assume that the observation
            # times are identical, if they are equal in number. Verifying the time array lengths
            # will avoid crashes due to interrupted simulations with incomplete data files.
            obs_data = (
                obs_0[0],
                (
                    np.asarray(obs_2[1]) - np.asarray(obs_0[1]) * np.asarray(obs_1[1])
                ).tolist(),
            )
            s_tex_label = (
                f"\\sigma^{s_obs_name[0]}_{{{q_indices[0]}}} "
                f"\\sigma^{s_obs_name[1]}_{{{q_indices[1]}}}"
            )
    return obs_data, s_tex_label


def prepare_2q_correlation_matrix(
    result: dict, s_obs_name: str, t: float, n_qubits: int
) -> (np.ndarray, str):
    """
    Prepare the data used for plotting the matrix of connected correlation values of one type for all.
    qubits. The connected correlation is defined as a 2Q observable from which the product of the two
    corresponding 1Q observables is subtracted.

    Args:
            result: A dictionary from which the observables are taken.
            s_obs_name: The name of the specific observable, used a key into the relevant observables
                    dict, and also in formatting the descriptive tex label of the data.
            t: The simulation time for which the data is to be calculated.
            n_qubits: The number of qubits in the simulation.

    Returns:
            A tuple with the following two entries:
                    obs_data: A matrix with the connected correlation function of all qubits at time `t`.
                    s_tex_label: A formatted tex label for the data - note that this string does not
                            indicate that the data is of 2Q connected correlation (the subtraction of
                            the 1Q product).
    """
    obs_1q_dict = result["obs-1q"]
    obs_2q_dict = result["obs-2q"]
    obs_data = np.full(shape=(n_qubits, n_qubits), dtype=float, fill_value=np.nan)
    s_obs_name = s_obs_name.lower()
    r_qubits = range(n_qubits)
    # Below we verify that all required 1Q and 2Q observables have complete data available.
    if obs_1q_dict is not None and obs_2q_dict is not None:
        for i in r_qubits:
            for j in r_qubits:
                if i == j:
                    continue
                obs_0 = obs_1q_dict.get((s_obs_name[0], (i,)))
                obs_1 = obs_1q_dict.get((s_obs_name[1], (j,)))
                obs_2 = obs_2q_dict.get((s_obs_name, (i, j)))
                if (
                    obs_0 is not None
                    and obs_1 is not None
                    and obs_2 is not None
                    and len(obs_0[0]) == len(obs_1[0])
                    and len(obs_0[0]) == len(obs_2[0])
                ):
                    # The data all comes from one simulation, so we can safely assume that the time
                    # arrays are identical, if they are equal in number. Verifying the time array lengths
                    # will avoid crashes due to interrupted simulations with incomplete data files.
                    try:
                        t_index = obs_0[0].index(t)
                        obs_data[i, j] = (
                            obs_2[1][t_index] - obs_0[1][t_index] * obs_1[1][t_index]
                        )
                    except ValueError:
                        pass
    s_tex_label = f"\\sigma^{s_obs_name[0]}_{{i}}\\sigma^{s_obs_name[1]}_{{j}}"
    return obs_data, s_tex_label


def plot_curves(
    obs_data_list: List[Tuple[Any, Any]],
    tex_labels: Optional[List[str]] = None,
    s_title="",
    ax=None,
    fontsize=16,
    line_styles: Optional[Sequence] = None,
    linewidth=3,
):
    """
    Plot multiple curves of simulation observables.

    Args:
            obs_data_list: A list of obs_data tuples of lists, from which the plotted data is taken.
            tex_labels: A list of descriptive tex label corresponding to the data.
            s_title: An optional plot title.
            ax: An optional axis object. If None, a new figure is created.
            fontsize: The fontsize to use in the figure.
            line_styles: A list with line styles iterated (periodically) for the curves. If None,
                    the default file-level member LINDBLADMPO_LINE_STYLES is used.
            linewidth: The plot line widths.

    Returns:
            An axis object (either the one passed as an argument, or a newly created one).
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))
    plt.rcParams.update({"font.size": fontsize})
    if line_styles is None:
        line_styles = LINDBLADMPO_LINE_STYLES
    n_styles = len(line_styles)
    for i_curve, obs_data in enumerate(obs_data_list):
        s_label = tex_labels[i_curve] if tex_labels is not None else None
        ax.plot(
            obs_data[0],
            obs_data[1],
            label=s_label,
            linestyle=line_styles[i_curve % n_styles],
            linewidth=linewidth,
        )
    if tex_labels is not None:
        ax.legend(fontsize=fontsize)
    ax.set_xlabel("$t$", fontsize=fontsize)
    if s_title != "":
        ax.set_title(s_title)
    return ax


def prepare_1q_space_time_data(
    parameters: dict,
    result: dict,
    s_obs_name: str,
    qubits: Optional[Union[List[int], np.ndarray]] = None,
    n_t_ticks=10,
    t_ticks_round=3,
    t_init: Optional[float] = None,
    t_final: Optional[float] = None,
) -> (np.ndarray, np.ndarray, np.ndarray, np.ndarray):
    """
    Prepare the data used for plotting a space-time (qubit-time) diagram of a single-qubit observable.

    Args:
            parameters: A dictionary from which the basic time parameters are taken.
            result: A dictionary from which the observables are taken.
            s_obs_name: The name of the specific single-qubit observable, used a key into the relevant
                    observables dict, and also in formatting the descriptive tex label of the data.
            qubits: An optional list of qubits to include in the data. If None, all qubits are used.
            n_t_ticks: The number of labeled major tick marks to generate for the time axis.
            t_ticks_round: The number of digits to round time axis tick labels to.
            t_init: An optional initial time - if None, it is taken from the parameters dict.
            t_final: An optional final time - if None, it is taken from the parameters dict.

    Returns:
            A tuple with the following four entries:
                    data: An array of the times for which the solver parameters indicate output is to
                            be evaluated for (based on the `t_init`, `t_final`, and `tau` parameters.
                    t_tick_indices: The indices of the time axis tick marks.
                    t_tick_labels: The formatted labels of the time axis tick marks.
                    qubits: The qubits used in the data.
    """

    t_eval, t_tick_indices, t_tick_labels, n_t_steps = prepare_time_data(
        parameters, n_t_ticks, t_ticks_round, t_init, t_final
    )
    if qubits is None:
        # Generate a default full 1Q matrix
        N = parameters["N"]
        qubits = np.arange(N)
    n_qubits = len(qubits)
    data = np.full(shape=(n_qubits, n_t_steps), dtype=float, fill_value=np.nan)
    for i_q, qubit in enumerate(qubits):
        obs_data, s_tex_label = prepare_curve_data(
            result, "obs-1q", s_obs_name, (qubit,)
        )
        if obs_data is not None:
            data[i_q, :] = obs_data[1][0:n_t_steps]
    return data, t_tick_indices, t_tick_labels, qubits


def prepare_2q_space_time_data(
    parameters: dict,
    result: dict,
    s_obs_name: str,
    qubit_0: Optional[int] = None,
    qubit_1: Optional[int] = None,
    n_t_ticks=10,
    t_ticks_round=3,
    t_init: Optional[float] = None,
    t_final: Optional[float] = None,
) -> (np.ndarray, np.ndarray, np.ndarray, np.ndarray):
    """
    Prepare the data used for plotting a space-time (qubit-time) diagram of a two-qubit observable.

    Args:
            parameters: A dictionary from which the basic time parameters are taken.
            result: A dictionary from which the observables are taken.
            s_obs_name: The name of the specific single-qubit observable, used a key into the relevant
                    observables dict, and also in formatting the descriptive tex label of the data.
            qubit_0: An optional index of the first qubit of the 2Q observable.
                    Exactly one of the arguments qubit_0 and qubit_1 must be None.
            qubit_1: An optional index of the second qubit of the 2Q observable.
                    Exactly one of the arguments qubit_0 and qubit_1 must be None.
            n_t_ticks: The number of labeled major tick marks to generate for the time axis.
            t_ticks_round: The number of digits to round time axis tick labels to.
            t_init: An optional initial time - if None, it is taken from the parameters dict.
            t_final: An optional final time - if None, it is taken from the parameters dict.

    Returns:
            A tuple with the following four entries:
                    data: An array with space-time (qubit-time) data.
                    t_tick_indices: The indices of the time axis tick marks.
                    t_tick_labels: The formatted labels of the time axis tick marks.
                    qubits: The qubits used in the data.
    """
    _, t_tick_indices, t_tick_labels, n_t_steps = prepare_time_data(
        parameters, n_t_ticks, t_ticks_round, t_init, t_final
    )
    N = parameters["N"]
    qubits = np.arange(N)
    b_q0 = qubit_0 is not None
    b_q1 = qubit_1 is not None
    if b_q0 and b_q1 or (not b_q0 and not b_q1):
        raise Exception(
            "Exactly one of the arguments qubit_0 and qubit_1 must be None."
        )
    n_qubits = len(qubits)
    data = np.full(shape=(n_qubits, n_t_steps), dtype=float, fill_value=np.nan)
    for i_q, qubit in enumerate(qubits):
        if not b_q0:
            qubits_pair = (qubit, qubit_1)
        else:  # else can be used according to the verification above
            qubits_pair = (qubit_0, qubit)
        obs_data, s_tex_label = prepare_curve_data(
            result, "obs-2q", s_obs_name, qubits_pair
        )
        if obs_data is not None:
            data[i_q, :] = obs_data[1][0:n_t_steps]
    return data, t_tick_indices, t_tick_labels, qubits


def prepare_2q_matrix_data(
    parameters: dict, result: dict, s_obs_name: str, t: Optional[float] = None
) -> (np.ndarray, np.ndarray):
    """
    Prepare the data used for plotting a two-qubit connected correlation matrix.

    Args:
            parameters: A dictionary from which the basic time parameters are taken.
            result: A dictionary from which the observables are taken.
            s_obs_name: The name of the specific two-qubit observable, used a key into the relevant
                    observables dict, and also in formatting the descriptive tex label of the data.
            t: An optional time for which to take the data. If unspecified the final time is used.

    Returns:
            A tuple with the following two entries:
                    data: An array with qubit-qubit data.
                    qubits: The qubits used in the data.
    """
    N = parameters["N"]
    if t is None:
        t = parameters["t_final"]
    qubits = np.arange(N)
    n_qubits = len(qubits)
    data, s_tex_label = prepare_2q_correlation_matrix(result, s_obs_name, t, n_qubits)
    return data, qubits


def plot_1q_space_time(
    data,
    s_obs_name: str,
    qubits,
    t_tick_indices,
    t_tick_labels,
    ax=None,
    fontsize=16,
    b_save_figures=True,
    s_file_prefix="",
):
    """
    Plot a single-qubit space-time (qubit-time) diagram.

    Args:
            data: The data to plot, prepared b a call to `prepare_1q_space_time_data()`.
            s_obs_name: The name of the single-qubit observable, used in formatting the descriptive
                    tex label of the data and the saved file name.
            qubits: The qubits used in the data.
            t_tick_indices: The indices of the time axis tick marks.
            t_tick_labels: The formatted labels of the time axis tick marks.
            ax: An optional axis object. If None, a new figure is created.
            fontsize: The fontsize to use in the figure.
            b_save_figures: Whether to save the plotted figure to file.
            s_file_prefix: An optional path and file name prefix for the saved figures.

    """
    s_obs_name = s_obs_name.lower()
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))
    plt.rcParams.update({"font.size": fontsize})
    im = ax.imshow(data, interpolation="none", aspect="auto")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    # plt.colorbar(im)
    ax.set_xlabel("$t$", fontsize=fontsize)
    ax.set_xticks(t_tick_indices)
    ax.set_xticklabels(t_tick_labels, fontsize=fontsize)
    ax.set_yticks(qubits)
    ax.set_yticklabels(qubits, fontsize=fontsize)
    ax.set_ylabel("$j$", fontsize=fontsize)
    s_title = f"$\\langle\\sigma^{s_obs_name}(t)\\rangle$"
    ax.set_title(s_title)
    s_file_label = f"sigma_{s_obs_name}.st"
    _save_fig(b_save_figures, s_file_prefix, s_file_label)


def plot_2q_correlation_matrix(
    data,
    s_obs_name: str,
    t: float,
    qubits,
    ax=None,
    fontsize=16,
    b_save_figures=True,
    s_file_prefix="",
    s_title=None,
):
    """
    Plot a two-qubit connected correlation matrix figure.

    Args:
            data: The data to plot, prepared b a call to `prepare_2q_matrix_data()`.
            s_obs_name: The name of the two-qubit observable, used in formatting the descriptive
                    tex label of the data and the saved file name.
            t: An optional time for which to take the data. If unspecified the final time is used.
            qubits: The qubits used in the data.
            ax: An optional axis object. If None, a new figure is created.
            fontsize: The fontsize to use in the figure.
            b_save_figures: Whether to save the plotted figure to file.
            s_file_prefix: An optional path and file name prefix for the saved figures.
            s_title: An optional title for the figure. If empty, a default title is formatted.
    """
    s_obs_name = s_obs_name.lower()
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 7))
    plt.rcParams.update({"font.size": fontsize})
    im = ax.imshow(data, interpolation="none", aspect="equal")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3.5%", pad=0.15)
    plt.colorbar(im, cax=cax)
    # plt.colorbar(im)
    n_qubits = len(qubits)
    if n_qubits <= 10:
        qubit_ticks = qubits
    else:
        step = 1 + (n_qubits // 10)
        qubit_ticks = qubits[0:n_qubits:step]
    ax.set_xticks(qubit_ticks)
    ax.set_xticklabels(qubit_ticks, fontsize=fontsize)
    ax.set_yticks(qubit_ticks)
    ax.set_yticklabels(qubit_ticks, fontsize=fontsize)
    ax.set_xlabel("$j$", fontsize=fontsize)
    ax.set_ylabel("$i$", fontsize=fontsize)
    s_tex_label = f"\\sigma^{s_obs_name[0]}_{{i}}\\sigma^{s_obs_name[1]}_{{j}}"
    if s_title is None:
        s_title = f"$\\langle{s_tex_label}(t={t})\\rangle_{{c}}$"
    ax.set_title(s_title)
    s_file_label = f"sigma_{s_obs_name[0]}.sigma_{s_obs_name[1]}.c.t={t}"
    _save_fig(b_save_figures, s_file_prefix, s_file_label)


def plot_full_1q_space_time(
    parameters: dict,
    result: dict,
    s_obs_name: str,
    ax=None,
    fontsize=16,
    b_save_figures=True,
    s_file_prefix="",
):
    """
    Prepare the data and plot a single-qubit space-time (qubit-time) diagram.

    Args:
            parameters: A dictionary from which the basic time parameters are taken.
            result: A dictionary from which the observables are taken.
            s_obs_name: The name of the single-qubit observable, used in formatting the descriptive
                    tex label of the data and the saved file name.
            ax: An optional axis object. If None, a new figure is created.
            fontsize: The fontsize to use in the figure.
            b_save_figures: Whether to save the plotted figure to file.
            s_file_prefix: An optional path and file name prefix for the saved figures.
    """
    data, t_tick_indices, t_tick_labels, qubits = prepare_1q_space_time_data(
        parameters, result, s_obs_name
    )
    plot_1q_space_time(
        data,
        s_obs_name,
        qubits,
        t_tick_indices,
        t_tick_labels,
        ax,
        fontsize,
        b_save_figures,
        s_file_prefix,
    )


def plot_full_2q_correlation_matrix(
    parameters: dict,
    result: dict,
    s_obs_name: str,
    t: Optional[float] = None,
    ax=None,
    fontsize=16,
    b_save_figures=True,
    s_file_prefix="",
    s_title = None,
):
    """
    Prepare the data and plot a two-qubit connected correlation matrix figure.

    Args:
            parameters: A dictionary from which the basic time parameters are taken.
            result: A dictionary from which the observables are taken.
            s_obs_name: The name of the two-qubit observable, used to obtain the data, and in formatting
                    the descriptive tex label of the data and the saved file name.
            t: An optional time for which to take the data. If unspecified the final time is used.
            ax: An optional axis object. If None, a new figure is created.
            fontsize: The fontsize to use in the figure.
            b_save_figures: Whether to save the plotted figure to file.
            s_file_prefix: An optional path and file name prefix for the saved figures.
            s_title: An optional title for the figure. If empty, a default title is formatted.
    """
    if t is None:
        t = parameters["t_final"]
    data, qubits = prepare_2q_matrix_data(parameters, result, s_obs_name, t)
    plot_2q_correlation_matrix(
        data,
        s_obs_name,
        t,
        qubits,
        ax,
        fontsize,
        b_save_figures,
        s_file_prefix,
        s_title,
    )


def plot_1q_obs_curves(
    parameters: dict,
    result: dict,
    s_obs_name: str,
    qubits: Optional[List[int]] = None,
    ax=None,
    fontsize=16,
    b_save_figures=True,
    s_file_prefix="",
    s_title=None,
):
    """
    Prepare the data and plot a single-qubit observable vs. time for multiple qubits.

    Args:
            parameters: A dictionary from which the basic time parameters are taken.
            result: A dictionary from which the observables are taken.
            s_obs_name: The name of the observable, used to obtain the data, and in formatting
                    the descriptive tex label of the data and the saved file name.
            qubits: The qubits to take for plotting.
            ax: An optional axis object. If None, a new figure is created.
            fontsize: The fontsize to use in the figure.
            b_save_figures: Whether to save the plotted figure to file.
            s_file_prefix: An optional path and file name prefix for the saved figures.
            s_title: An optional title for the figure. If empty, a default title is formatted.
    """
    obs_data_list = []
    tex_labels = []
    s_obs_name = s_obs_name.lower()
    for qubit in qubits:
        obs_data, s_tex_label = prepare_curve_data(
            result, "obs-1q", s_obs_name, (qubit,)
        )
        obs_data_list.append(obs_data)
        tex_labels.append(f"$\\langle{s_tex_label}(t)\\rangle$")
    if s_title is None:
        s_title = f"$\\langle\\sigma^{s_obs_name}_j(t)\\rangle$"
    ax = plot_curves(obs_data_list, tex_labels, s_title, ax, fontsize)
    _, t_tick_indices, t_tick_labels, _ = prepare_time_data(parameters)
    # ax.set_xticks(t_tick_indices)
    ax.set_xticks(t_tick_labels)
    ax.set_xticklabels(t_tick_labels, fontsize=fontsize)
    s_file_label = f"sigma_{s_obs_name}"
    _save_fig(b_save_figures, s_file_prefix, s_file_label)


def plot_2q_obs_curves(
    parameters: dict,
    result: dict,
    s_obs_name: str,
    qubit_pairs: Optional[List[Tuple[int, int]]] = None,
    ax=None,
    fontsize=16,
    b_save_figures=True,
    s_file_prefix="",
    s_title=None,
):
    """
    Prepare the data and plot a two-qubit observable vs. time for multiple qubits.

    Args:
            parameters: A dictionary from which the basic time parameters are taken.
            result: A dictionary from which the observables are taken.
            s_obs_name: The name of the observable, used to obtain the data, and in formatting
                    the descriptive tex label of the data and the saved file name.
            qubit_pairs: The qubit pairs to take for plotting.
            ax: An optional axis object. If None, a new figure is created.
            fontsize: The fontsize to use in the figure.
            b_save_figures: Whether to save the plotted figure to file.
            s_file_prefix: An optional path and file name prefix for the saved figures.
            s_title: An optional title for the figure. If empty, a default title is formatted.
    """
    obs_data_list = []
    tex_labels = []
    s_obs_name = s_obs_name.lower()
    for q_pair in qubit_pairs:
        obs_data, s_tex_label = prepare_curve_data(result, "obs-2q", s_obs_name, q_pair)
        if obs_data is not None:
            obs_data_list.append(obs_data)
            tex_labels.append(f"$\\langle{s_tex_label}(t)\\rangle$")
    if s_title is None:
        s_title = (
            f"$\\langle\\sigma^{s_obs_name[0]}_i\\sigma^{s_obs_name[1]}_j(t)\\rangle$"
        )
    ax = plot_curves(obs_data_list, tex_labels, s_title, ax, fontsize)
    _, t_tick_indices, t_tick_labels, _ = prepare_time_data(parameters)
    # ax.set_xticks(t_tick_indices)
    ax.set_xticks(t_tick_labels)
    ax.set_xticklabels(t_tick_labels, fontsize=fontsize)
    s_file_label = f"sigma_{s_obs_name[0]}.sigma_{s_obs_name[1]}"
    _save_fig(b_save_figures, s_file_prefix, s_file_label)


def plot_2q_correlation_curves(
    parameters: dict,
    result: dict,
    s_obs_name: str,
    qubit_pairs: Optional[List[Tuple[int, int]]] = None,
    ax=None,
    fontsize=16,
    b_save_figures=True,
    s_file_prefix="",
    s_title=None,
):
    """
    Prepare the data and plot a two-qubit connected correlation curve vs. time for multiple qubits.

    Args:
            parameters: A dictionary from which the basic time parameters are taken.
            result: A dictionary from which the observables are taken.
            s_obs_name: The name of the observable, used to obtain the data, and in formatting
                    the descriptive tex label of the data and the saved file name.
            qubit_pairs: The qubit pairs to take for plotting.
            ax: An optional axis object. If None, a new figure is created.
            fontsize: The fontsize to use in the figure.
            b_save_figures: Whether to save the plotted figure to file.
            s_file_prefix: An optional path and file name prefix for the saved figures.
            s_title: An optional title for the figure. If empty, a default title is formatted.
    """
    obs_data_list = []
    tex_labels = []
    s_obs_name = s_obs_name.lower()
    for q_pair in qubit_pairs:
        obs_data, s_tex_label = prepare_2q_correlation_data(result, s_obs_name, q_pair)
        if obs_data is not None:
            obs_data_list.append(obs_data)
            tex_labels.append(f"$\\langle{s_tex_label}(t)\\rangle_{{c}}$")
    if s_title is None:
        s_title = f"$\\langle\\sigma^{s_obs_name[0]}_i\\sigma^{s_obs_name[1]}_j(t)\\rangle_{{c}}$"
    ax = plot_curves(obs_data_list, tex_labels, s_title, ax, fontsize)
    _, t_tick_indices, t_tick_labels, _ = prepare_time_data(parameters)
    # ax.set_xticks(t_tick_indices)
    ax.set_xticks(t_tick_labels)
    ax.set_xticklabels(t_tick_labels, fontsize=fontsize)
    s_file_label = f"sigma_{s_obs_name[0]}.sigma_{s_obs_name[1]}.c"
    _save_fig(b_save_figures, s_file_prefix, s_file_label)


def plot_global_obs_curve(
    parameters: dict,
    result: dict,
    s_obs_name: str,
    ax=None,
    fontsize=16,
    b_save_figures=True,
    s_file_prefix="",
):
    """
    Prepare the data and plot a global quantity curve vs. time.

    Args:
            parameters: A dictionary from which the basic time parameters are taken.
            result: A dictionary from which the observables are taken.
            s_obs_name: The name of the observable, used to obtain the data, and in formatting
                    the descriptive tex label of the data and the saved file name.
            ax: An optional axis object. If None, a new figure is created.
            fontsize: The fontsize to use in the figure.
            b_save_figures: Whether to save the plotted figure to file.
            s_file_prefix: An optional path and file name prefix for the saved figures.
    """
    s_obs_name = s_obs_name.lower()
    obs_data, s_tex_label = prepare_curve_data(result, "global", s_obs_name, ())
    s_title = f"${s_tex_label}(t)$"
    ax = plot_curves([obs_data], [s_title], s_title, ax, fontsize)
    _, t_tick_indices, t_tick_labels, _ = prepare_time_data(parameters)
    # ax.set_xticks(t_tick_indices)
    ax.set_xticks(t_tick_labels)
    ax.set_xticklabels(t_tick_labels, fontsize=fontsize)
    s_file_label = s_obs_name
    _save_fig(b_save_figures, s_file_prefix, s_file_label)


def _save_fig(b_save_figures, s_file_prefix, s_file_label):
    if b_save_figures:
        if s_file_prefix != "":
            s_file_label = "." + s_file_label
        plt.savefig(s_file_prefix + s_file_label + ".png")
