#include <Siv3D.hpp> // OpenSiv3D v0.6.0

#include "include/GraphDrawing.hpp"

struct Station
{
	String name;
	Vec2 pos;

	Station() = default;

	Station(const String& name, const Vec2& pos)
		: name(name), pos(pos)
	{}
};

struct EdgeInfo
{
	Array<String> lineNames;
};

class RailwayVisualizer : public IGraphVisualizer
{
public:

	RailwayVisualizer() = default;

	virtual ~RailwayVisualizer() = default;

	void drawNode(const Vec2& pos, GraphEdge::IndexType nodeIndex) const override
	{
		const auto text = m_font(m_stationsNames.find(nodeIndex)->second);
		const auto regionRect = text.regionAt(m_font.fontSize(), pos).stretched(m_font.fontSize(), 0);

		const auto rectRounded = regionRect.rounded(regionRect.h * 0.5);

		rectRounded.draw(Palette::White);
		rectRounded.drawFrame(2, Palette::Black);

		text.drawAt(pos, Palette::Black);
	}

	void drawEdge(const Line& line, GraphEdge::IndexType index0, GraphEdge::IndexType index1) const override
	{
		const Point edgeKey{ Min(index0, index1), Max(index0, index1) };

		const auto& edgeInfo = m_edgeInfos.find(edgeKey)->second;

		const Vec2 v = line.vector().normalize();
		const Vec2 n = Vec2(-v.y, v.x);

		const double interval = 7;

		for (const auto& [i, lineName] : Indexed(edgeInfo.lineNames))
		{
			const auto color = m_lineColors.find(lineName)->second;

			line.movedBy(n * interval * static_cast<int>(i)).draw(5, color);
		}
	}

	Font m_font = Font{ 12 };

	HashTable<String, Color> m_lineColors;

	HashTable<GraphEdge::IndexType, String> m_stationsNames;

	HashTable<Point, EdgeInfo> m_edgeInfos;
};

void Main()
{
	TOMLReader reader(U"./trainData.toml");

	RailwayVisualizer visualizer;

	visualizer.m_lineColors = {
		{U"JR山手線", Color(139, 198, 62)},
		{U"JR中央・総武線", Color(254, 211, 4)},
		{U"JR埼京線", Color(0, 133, 64)},
		{U"JR京浜東北線", Color(0, 174, 239)},
		{U"JR湘南新宿ライン", Color(237, 108, 0)},
		{U"宇都宮線", Color(240, 134, 0)},
	};

	const Rect clipRect = Rect(1500, 2000).setCenter(0, 0);

	HashTable<int, Vec2> initialNodePositions;

	for (const auto& stationData : reader[U"Station"].tableArrayView())
	{
		const auto id = stationData[U"id"].get<int32>();
		const auto name = stationData[U"name"].getString();
		const auto pos = stationData[U"pos"].get<Vec2>();

		if (clipRect.contains(pos))
		{
			initialNodePositions[id] = pos;
			visualizer.m_stationsNames[id] = name;
		}
	}

	Array<GraphEdge> graphEdges;

	for (const auto& lineData : reader[U"Line"].tableArrayView())
	{
		const auto name = lineData[U"name"].getString();

		for (const auto& edgeData : lineData[U"Edge"].tableArrayView())
		{
			const auto idA = edgeData[U"idA"].get<int32>();
			const auto idB = edgeData[U"idB"].get<int32>();

			if (initialNodePositions.find(idA) != initialNodePositions.end() && initialNodePositions.find(idB) != initialNodePositions.end())
			{
				const Point edgeKey{ Min(idA, idB), Max(idA, idB) };

				graphEdges.push_back({ edgeKey.x, edgeKey.y });

				visualizer.m_edgeInfos[edgeKey].lineNames.push_back(name);
			}
		}
	}

	GraphSet graph(graphEdges);

	const ForceDirectedConfig config{ .K = 30 };

	LayoutForceDirected layout(graph[0], config, initialNodePositions);

	Scene::SetBackground(Palette::White);

	Window::Resize(800, 800);

	while (System::Update())
	{
		layout.update(1);

		layout.setDrawArea(Rect(800, 800));

		layout.draw(visualizer);
	}
}
