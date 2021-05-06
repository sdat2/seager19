"""ani.py - A set of functions to animate particular results.

animate_xr_da - animation for individual xr.DataArray.

animate_prediction - plots inputs and outputs of old model.

"""
import os
from typing import Callable
import numpy as np
import xarray as xr
from tqdm import tqdm
import matplotlib
import matplotlib.pyplot as plt
import imageio
from src.plot_settings import label_subplots, ps_defaults, time_title, cmap
from src.utils import timeit, fix_calendar


def _rdict(index: int) -> dict:
    """Returns renaming dict for xarray.DataArray.

    Made to reformat the output datarrays of the Fortran
    ocean model used.

    Args:
        index (int): index on coords

    Returns:
        dict: renaming dict.

    """
    return {
        "T_0" + str(index): "time",
        "Y_0" + str(index): "y",
        "X_0" + str(index): "x",
        "L_0" + str(index): "Z",
    }


@timeit
def animate_ds(ds: xr.Dataset, file_name: str, output_dir: str) -> None:
    """Animate the `xarray.Dataset`.

    Args:
        ds (xr.Dataset): xarray.Dataset to animate the variables of.
        file_name (str): Name of dataset to be associated with the animations.
        output_dir (str): Full path to output directory to put the animations in.

    """
    ps_defaults(use_tex=False, dpi=200)
    cmap_d = {
        "DYN_PRES": "delta",
        "SST_QFLX": "delta",
        "SST_SST": "sst",
        "qflx": "delta",
        "QFLX": "delta",
        "TDEEP_HTHERM": "sst",
        "TDEEP_TDEEP": "sst",
        "TDEEP_HMODEL": "sst",
    }
    unit_d = {"SST_SST": r"$^{\circ}$C"}
    for y in ds.variables:
        y = str(y)
        if "X_" not in y:
            if "Y_" not in y:
                if "L_" not in y:
                    if "T_" not in y or "SST" in y:
                        if "GRID" != y:
                            print(y)
                            da = ds[y]
                            if (
                                "T_01" in da.coords
                                or "T_02" in da.coords
                                or "T_03" in da.coords
                                or "T_04" in da.coords
                            ):
                                for key in da.coords:
                                    num = str(key)[3]
                                da = da.rename(_rdict(num))
                            if y in unit_d:
                                da.attrs["units"] = unit_d[y]
                            da.coords["x"].attrs["units"] = r"$^{\circ}$E"
                            da.coords["y"].attrs["units"] = r"$^{\circ}$N"
                            da = da.where(da != 0.0).isel(Z=0)
                            da = fix_calendar(da, timevar="time")
                            if "variable" in da.dims:
                                da = da.isel(variable=0)
                            da = da.rename(y)
                            if y in unit_d:
                                da.attrs["units"] = unit_d[y]
                            da.attrs["long_name"] = y
                            da.attrs["name"] = y
                            animate_xr_da(
                                da,
                                video_path=os.path.join(
                                    output_dir, file_name + "_" + y + ".gif"
                                ),
                                vcmap=cmap_d[y],
                            )


@timeit
def animate_xr_da(
    xr_da: xr.DataArray,
    video_path: str = "output.mp4",
    vcmap: any = cmap("sst"),
) -> None:
    """Animate an `xr.DataArray`.

    Args:
        xr_da (xr.DataArray): input xr.DataArray.
        video_path (str, optional): Video path. Defaults to "output.mp4".
        vcmap (any, optional): cmap for variable. Defaults to cmap("sst").

    """
    ps_defaults(use_tex=False, dpi=200)
    balanced_colormap = False

    if isinstance(vcmap, str):
        if vcmap == "delta":
            balanced_colormap = True
        vcmap = cmap(vcmap)

    assert isinstance(vcmap, matplotlib.colors.LinearSegmentedColormap)

    def gen_frame_func(
        xr_da: xr.DataArray,
    ) -> Callable:
        """Create imageio frame function for `xarray.DataArray` visualisation.

        Args:
            x_da (xr.DataArray): input xr.DataArray.

        Returns:
            make_frame (Callable): function to create each frame.

        """
        vmin = xr_da.min(skipna=True)
        vmax = xr_da.max(skipna=True)
        if balanced_colormap:
            vmin, vmax = [np.min([vmin, -vmax]), np.max([vmax, -vmin])]

        def make_frame(index: int) -> np.array:
            """Make an individual frame of the animation.

            Args:
                index (int): The time index.

            Returns:
                image (np.array): np.frombuffer output that can be fed into imageio

            """
            fig, ax1 = plt.subplots(1, 1)

            xr_da.isel(time=index).plot.imshow(ax=ax1, cmap=vcmap, vmin=vmin, vmax=vmax)
            time_title(ax1, xr_da.time.values[index])
            plt.tight_layout()

            fig.canvas.draw()
            image = np.frombuffer(fig.canvas.tostring_rgb(), dtype="uint8")
            image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
            plt.close()

            return image

        return make_frame

    def xarray_to_video(
        xr_da: xr.DataArray,
        video_path: str,
        fps: int = 5,
    ) -> None:
        """Generate video of an `xarray.DataArray`.

        Args:
            xr_da (xr.DataArray): input xarray.DataArray
            video_path (str, optional): output path to save.
            fps (int, optional): frames per second.

        """
        video_indices = list(range(len(xr_da.time.values)))
        make_frame = gen_frame_func(xr_da)
        imageio.mimsave(
            video_path,
            [make_frame(index) for index in tqdm(video_indices, desc=video_path)],
            fps=fps,
        )
        print("Video " + video_path + " made.")

    xarray_to_video(xr_da, video_path, fps=5)


@timeit
def animate_prediction(
    x_da: xr.DataArray,
    y_da: xr.DataArray,
    pred_da: xr.DataArray,
    video_path: str = "joint_val.mp4",
) -> None:
    """This function animates the inputs, labels, predictions.

    Args:
        x_da (xr.DataArray): 3 or 6 bands, 4 seasons, 20 years
        y_da (xr.DataArray): 1 band, 20 years
        pred_da (xr.DataArray): 1 band, 20 years
        video_path (str, optional): relative text path to output mp4 file.
            Defaults to "joint_val.mp4".

    """

    def gen_frame_func(
        x_da: xr.DataArray, y_da: xr.DataArray, pred_da: xr.DataArray
    ) -> Callable:
        """Create imageio frame function for xarray.DataArray visualisation.

        Args:
            x_da (xr.DataArray): 3 or 6 bands, 4 seasons, 20 years
            y_da (xr.DataArray): 1 band, 20 years
            pred_da (xr.DataArray): 1 band, 20 years

        Returns:
            make_frame (Callable): function to create each frame.

        """

        def make_frame(index: int) -> np.array:
            """Make an individual frame of the animation.

            Args:
                index (int): The time index.

            Returns:
                image (np.array): np.frombuffer output that can be fed into imageio

            """
            if len(x_da.band.values) == 3:
                fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(
                    3, 2, figsize=(10, 10)
                )
            elif len(x_da.band.values) == 6:
                fig, (
                    (ax1, ax2),
                    (ax1b, ax2b),
                    (ax3, ax4),
                    (ax3b, ax4b),
                    (ax5, ax6),
                ) = plt.subplots(5, 2, figsize=(10, 17))
            else:
                assert False

            x_da.isel(year=index, mn=0, band=slice(0, 3)).plot.imshow(ax=ax1)
            x_da.isel(year=index, mn=1, band=slice(0, 3)).plot.imshow(ax=ax2)
            x_da.isel(year=index, mn=2, band=slice(0, 3)).plot.imshow(ax=ax3)
            x_da.isel(year=index, mn=3, band=slice(0, 3)).plot.imshow(ax=ax4)

            for ax in [ax1, ax2, ax3, ax4]:
                ax.set_xlabel("")

            if len(x_da.band.values) == 6:
                x_da.isel(year=index, mn=0, band=slice(3, 6)).plot.imshow(ax=ax1b)
                x_da.isel(year=index, mn=1, band=slice(3, 6)).plot.imshow(ax=ax2b)
                x_da.isel(year=index, mn=2, band=slice(3, 6)).plot.imshow(ax=ax3b)
                x_da.isel(year=index, mn=3, band=slice(3, 6)).plot.imshow(ax=ax4b)
                for ax in [ax1b, ax2b, ax3b, ax4b]:
                    ax.set_xlabel("")
                label_subplots(
                    [ax1, ax2, ax1b, ax2b, ax3, ax4, ax3b, ax4b, ax5, ax6], y_pos=1.07
                )
            else:
                label_subplots([ax1, ax2, ax3, ax4, ax5, ax6], y_pos=1.07)

            xr.DataArray(
                data=y_da.isel(year=index).values,
                dims=["y", "x"],
                coords=dict(
                    y=y_da.coords["y"].values,
                    x=y_da.coords["x"].values,
                    Y="esa_cci",
                ),
            ).plot(ax=ax5)

            xr.DataArray(
                data=pred_da.isel(year=index),
                dims=["y", "x"],
                coords=dict(
                    y=y_da.coords["y"].values,
                    x=y_da.coords["x"].values,
                    Y="predicted_classes",
                ),
            ).plot(ax=ax6)
            fig.canvas.draw()
            image = np.frombuffer(fig.canvas.tostring_rgb(), dtype="uint8")
            image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
            plt.close()
            return image

        return make_frame

    def xarray_to_video(
        x_da: xr.DataArray,
        y_da: xr.DataArray,
        pred_da: xr.DataArray,
        video_path: str,
        fps: int = 5,
    ) -> None:
        """Generate video of an `xarray.DataArray`.

        The full set of time coordinates of the datasets are used.

        Args:
            x_da (xr.DataArray): 3 or 6 bands, 4 seasons, 20 years
            y_da (xr.DataArray): 1 band, 20 years
            pred_da (xr.DataArray): 1 band, 20 years
            video_path (str, optional): relative text path to output mp4 file.
            fps (int, optional): frames per second.

        """
        video_indices = list(range(len(y_da.year.values)))
        make_frame = gen_frame_func(x_da, y_da, pred_da)
        imageio.mimsave(
            video_path,
            [make_frame(index) for index in tqdm(video_indices, desc=video_path)],
            fps=fps,
        )
        print("Video " + video_path + " made.")

    xarray_to_video(x_da, y_da, pred_da, video_path, fps=5)
