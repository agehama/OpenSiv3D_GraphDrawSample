#include <Siv3D.hpp> // OpenSiv3D v0.6

#include "include/GraphDrawing.hpp"

class TreeGraphData
{
public:

	TreeGraphData() = default;

	TreeGraphData(const JSON& json, const String& rootLabel, Optional<int> maxDepth = none, Optional<int> maxArrayWidth = none, Optional<int> maxObjectWidth = none)
	{
		loadJSON(json, rootLabel, maxDepth, maxArrayWidth, maxObjectWidth);
	}

	GraphEdge::IndexType newNode(const String& label, JSONValueType type, int depth)
	{
		const GraphEdge::IndexType index = static_cast<GraphEdge::IndexType>(m_nodes.size());
		m_nodes.push_back(label);
		m_types.push_back(type);
		m_depths.push_back(depth);

		return index;
	}

	void loadJSON(const JSON& json, const String& rootLabel, Optional<int> maxDepth = none, Optional<int> maxArrayWidth = none, Optional<int> maxObjectWidth = none)
	{
		m_nodes.clear();
		m_types.clear();
		m_depths.clear();
		m_edges.clear();

		loadJSONObject(json, newNode(rootLabel, JSONValueType::Object, 0), 1, maxDepth, maxArrayWidth, maxObjectWidth);
	}

	int maxDepth()const
	{
		return *std::max_element(m_depths.begin(), m_depths.end());
	}

	const Array<String>& nodes() const { return m_nodes; }

	const Array<JSONValueType>& types() const { return m_types; }

	const Array<int>& depths() const { return m_depths; }

	const Array<GraphEdge>& edges() const { return m_edges; }

private:

	void loadJSONPrimitive(const JSON& json, GraphEdge::IndexType parent, int depth)
	{
		String label;

		if (json.isString())
		{
			label = json.get<String>();
		}
		else if (json.isInteger())
		{
			label = Format(json.get<int>());
		}
		else if (json.isFloat())
		{
			label = Format(json.get<float>());
		}
		else if (json.isBool())
		{
			label = Format(json.get<bool>());
		}
		else if (json.isNull())
		{
			label = U"null";
		}
		else
		{
			label = U"empty";
		}

		const GraphEdge::IndexType nodeIndex = newNode(label, json.getType(), depth);

		m_edges.emplace_back(parent, nodeIndex);
	}

	void loadJSONArray(const JSON& json, GraphEdge::IndexType parentIndex, int depth, Optional<int> maxDepth, Optional<int> maxArrayWidth, Optional<int> maxObjectWidth)
	{
		if (maxDepth && maxDepth.value() <= depth)
		{
			m_edges.emplace_back(parentIndex, newNode(U"...", JSONValueType::Object, depth));

			return;
		}

		GraphEdge::IndexType i = 0;

		for (const auto& val : json.arrayView())
		{
			if (maxArrayWidth && maxArrayWidth.value() <= i)
			{
				m_edges.emplace_back(parentIndex, newNode(U"...", JSONValueType::Object, depth));

				return;
			}

			const GraphEdge::IndexType nodeIndex = newNode(Format(U"[", i, U"]"), val.getType(), depth);

			m_edges.emplace_back(parentIndex, nodeIndex);

			if (val.isArray())
			{
				loadJSONArray(val, nodeIndex, depth + 1, maxDepth, maxArrayWidth, maxObjectWidth);
			}
			else if (val.isObject())
			{
				loadJSONObject(val, nodeIndex, depth + 1, maxDepth, maxArrayWidth, maxObjectWidth);
			}
			else
			{
				loadJSONPrimitive(val, nodeIndex, depth + 1);
			}

			++i;
		}
	}

	void loadJSONObject(const JSON& json, GraphEdge::IndexType parentIndex, int depth, Optional<int> maxDepth, Optional<int> maxArrayWidth, Optional<int> maxObjectWidth)
	{
		if (maxDepth && maxDepth.value() <= depth)
		{
			m_edges.emplace_back(parentIndex, newNode(U"...", JSONValueType::Object, depth));

			return;
		}

		GraphEdge::IndexType i = 0;

		for (auto& [key, val] : json)
		{
			if (maxObjectWidth && maxObjectWidth.value() <= i)
			{
				m_edges.emplace_back(parentIndex, newNode(U"...", JSONValueType::Object, depth));

				return;
			}

			const GraphEdge::IndexType nodeIndex = newNode(key, JSONValueType::Object, depth);

			m_edges.emplace_back(parentIndex, nodeIndex);

			if (val.isArray())
			{
				loadJSONArray(val, nodeIndex, depth + 1, maxDepth, maxArrayWidth, maxObjectWidth);
			}
			else if (val.isObject())
			{
				loadJSONObject(val, nodeIndex, depth + 1, maxDepth, maxArrayWidth, maxObjectWidth);
			}
			else
			{
				loadJSONPrimitive(val, nodeIndex, depth + 1);
			}

			++i;
		}
	}

	Array<String> m_nodes;

	Array<JSONValueType> m_types;

	Array<int> m_depths;

	Array<GraphEdge> m_edges;
};

class TreeGraphVisualizer : public BasicGraphVisualizer
{
public:

	explicit TreeGraphVisualizer(TreeGraphData& graphData, const Font& font, const ColorF& fontColor, const HSV& rootNodeColor, const ColorF& edgeColor, double edgeThickness)
		: BasicGraphVisualizer(0, edgeThickness, ColorF{}, edgeColor)
		, m_font(font)
		, m_fontColor(fontColor)
		, m_rootNodeColor(rootNodeColor)
		, graphDataRef(graphData)
	{}

	void drawNode(const Vec2& pos, GraphEdge::IndexType index) const override
	{
		const double absScale = 1.0 / Graphics2D::GetMaxScaling();

		const auto& graphData = graphDataRef.get();

		const double fontSize = m_font.fontSize() * absScale;

		const auto resionRect = m_font(graphData.nodes()[index]).regionAt(fontSize, pos).stretched(fontSize, 0);

		const auto shadowOffset = Circular{ fontSize * 0.5, 150_deg }.fastToVec2();

		const HSV nodeColor{ m_rootNodeColor.a + graphData.depths()[index] * 30, m_rootNodeColor.s, m_rootNodeColor.v };

		resionRect.rounded(resionRect.h * 0.5).drawShadow(shadowOffset, fontSize * 0.5, 1.0 * absScale, ColorF(0.1, 0.4)).draw(nodeColor);

		// 文字の影
		m_font(graphData.nodes()[index]).drawAt(fontSize, pos + Vec2(1, 1) * absScale, ColorF(0.3));

		m_font(graphData.nodes()[index]).drawAt(fontSize, pos, m_fontColor);
	}

	void drawEdge(const Line& line, GraphEdge::IndexType index0, GraphEdge::IndexType index1) const override
	{
		const auto& graphData = graphDataRef.get();

		const double depthRate = (graphData.depths()[index0] + graphData.depths()[index1]) * 0.5 / graphData.maxDepth();

		const double thickness = Math::Lerp(m_edgeThickness, 1.0, depthRate) / Graphics2D::GetMaxScaling();

		line.draw(LineStyle::RoundCap, thickness, m_edgeColor);
	}

	Font m_font;

	Color m_fontColor;

	HSV m_rootNodeColor;

	std::reference_wrapper<TreeGraphData> graphDataRef;
};


void Main()
{
	// 読み込み上限の設定
	// JSONが大きすぎると時間がかかり見た目も複雑になるので上限をかける
	// 上限を超えて省略されたツリーは"..."で表示される
	const Optional<int> maxDepth = none; // 木の深さの上限
	const Optional<int> maxArrayWidth = 5; // 各Arrayの要素数の上限
	const Optional<int> maxObjectWidth = 10; // 各Objectの要素数の上限

	const Color backgroundColor{ U"#d4cdb8" };
	const Color edgeColor{ U"#595154" };
	const Color fontColor{ U"#e4d7bf" };
	const HSV rootNodeColor{ 0, 0.3, 0.7 };

	Window::Resize(1000, 1000);

	TreeGraphData graphData{ JSON::Load(U"example/json/config.json"), U"config.json" };

	const ForceDirectedConfig config{ .repulsiveExponent = 3 };

	LayoutForceDirected layout1{ graphData.edges(), config };

	Scene::SetBackground(backgroundColor);

	TreeGraphVisualizer visualizer{ graphData, Font{ 14, Typeface::Black }, fontColor, rootNodeColor, edgeColor, 15 };

	Camera2D camera{ Scene::Rect().center() };

	Window::SetTitle(U"JSONViewer -  [ドラッグドロップ]:JSONファイルを開く  [マウス]:カメラ操作  [Z]:文字サイズ小  [X]:文字サイズ大  [C]:カメラリセット  [Space]:再配置");

	while (System::Update())
	{
		if (DragDrop::HasNewFilePaths())
		{
			const auto paths = DragDrop::GetDroppedFilePaths();
			const auto& path = paths[0].path;

			if (const auto json = JSON::Load(path))
			{
				graphData.loadJSON(json, FileSystem::FileName(path), maxDepth, maxArrayWidth, maxObjectWidth);

				layout1.init(graphData.edges(), GetDefaultRNG(), config);

				camera = Camera2D{ Scene::Rect().center() };
			}
		}

		if (KeySpace.down())
		{
			layout1.reset();
		}

		if (KeySpace.down() || KeyC.down())
		{
			camera = Camera2D{ Scene::Rect().center() };
		}

		if (KeyZ.down())
		{
			visualizer.m_font = Font{ visualizer.m_font.fontSize() - 1, Typeface::Black };
		}

		if (KeyX.down())
		{
			visualizer.m_font = Font{ visualizer.m_font.fontSize() + 1, Typeface::Black };
		}

		{
			auto t = camera.createTransformer();

			layout1.update(16);

			layout1.setDrawArea(Scene::Rect().stretched(-50));

			layout1.draw(visualizer);
		}

		camera.update();

		camera.draw(Palette::White);
	}
}
