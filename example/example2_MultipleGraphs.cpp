#include <Siv3D.hpp> // OpenSiv3D v0.6

#include "include/GraphDrawing.hpp"

class LabelGraphVisualizer : public BasicGraphVisualizer
{
public:

	explicit LabelGraphVisualizer(const Font& font, ColorF fontColor, double nodeRadius = 10.0, double edgeThickness = 1.0, ColorF nodeColor = Palette::White, ColorF edgeColor = ColorF(0.8))
		: BasicGraphVisualizer{ nodeRadius, edgeThickness, nodeColor, edgeColor }
		, m_labelFont(font)
		, m_labelColor(fontColor)
	{}

	virtual ~LabelGraphVisualizer() = default;

	virtual void drawNode(const Vec2& pos, GraphEdge::IndexType nodeIndex) const override
	{
		pos.asCircle(m_nodeRadius).draw(m_nodeColor);
		m_labelFont(nodeIndex).drawAt(pos, m_labelColor);
	}

	Font m_labelFont;

	ColorF m_labelColor;
};

void Main()
{
	Window::Resize({ 640,640 });

	const GraphLoader loader{ U"./primitives.mtx" };

	const ForceDirectedConfig config{
		.startImmediately = StartImmediately::Yes,
		.tol = 0.001,
	};

	Reseed(0);
	auto layouts = Array<LayoutForceDirected>::IndexedGenerate(loader.size(), [&](size_t i) {
		return LayoutForceDirected{ loader[i], config };
		}
	);

	const int margin = 30;
	const int columns = 3;
	const Point rectSize = (Scene::Size() - Point{ margin, margin }*2) / columns;

	const Color backgroundColor{ Palette::White };
	const Color foregroundColor{ Palette::Black };
	const Color fontColor{ backgroundColor };

	Scene::SetBackground(backgroundColor);
	LabelGraphVisualizer visualizer{ Font{16, Typeface::Black }, fontColor, 12, 6, foregroundColor, foregroundColor };

	while (System::Update())
	{
		for (auto [i, layout] : IndexedRef(layouts))
		{
			const auto rect = Rect{ rectSize }.setPos(Point{ margin, margin } + rectSize * Point{ i % columns, i / columns });

			layout.setDrawArea(rect.stretched(-50));

			layout.draw(visualizer);

			rect.stretched(-25).drawFrame(5.0, foregroundColor);
		}
	}
}
