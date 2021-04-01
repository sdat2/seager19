"""ani.py - A set of functions to animate particular results.

animate_prediction - specifically designed to plot inputs and results
    from the xgboost algorithm over Chenobyl.

"""
from typing import Callable
import numpy as np
import xarray as xr
from tqdm import tqdm
import matplotlib.pyplot as plt
import imageio
from src.plot_settings import label_subplots, ps_defaults
from src.utils import timeit

ps_defaults(use_tex=False, dpi=200)


@timeit
def animate_prediction(
    x_da: xr.DataArray,
    y_da: xr.DataArray,
    pred_da: xr.DataArray,
    video_path: str = "joint_val.mp4",
) -> None:
    """This function animates the inputs, labels, predictions.

    Based on code originally from Tom Anderson: tomand@bas.ac.uk.

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
        """Generate video of an xarray.DataArray.
        The full set of time coordinates of the datasets are used.
        Args:
            x_da (xr.DataArray): 3 or 6 bands, 4 seasons, 20 years
            y_da (xr.DataArray): 1 band, 20 years
            pred_da (xr.DataArray): 1 band, 20 years
            video_path (str, optional): relative text path to output mp4 file.
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
