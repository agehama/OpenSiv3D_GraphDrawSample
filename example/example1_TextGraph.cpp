#include <Siv3D.hpp> // OpenSiv3D v0.6

#include "include/GraphDrawing.hpp"

void Main()
{
	Window::Resize(1280, 720);

	const GraphLoader loader{ U"./textGraph.mtx" };

	Reseed(0);
	LayoutForceDirected layout{ loader[0], ForceDirectedConfig{.repulsiveExponent = 3.5 } };

	const BasicGraphVisualizer visualizer{ 1.0, 0.0 };

	Camera2D camera{ Scene::Center() };

	const double rotate = -200_deg;
	const auto rotateMat = Mat3x2::Rotate(rotate, Scene::Center());

	const auto drawAreaRect = Scene::Rect().stretched(-50);

	while (System::Update())
	{
		{
			auto transformer = camera.createTransformer();

			Transformer2D transformer2{ rotateMat };

			layout.update(16);

			layout.setDrawArea(drawAreaRect, rotateMat);

			layout.draw(visualizer);
		}

		camera.update();

		camera.draw(Palette::Orange);
	}
}
