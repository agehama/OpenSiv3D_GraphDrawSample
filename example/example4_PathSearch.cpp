#include <Siv3D.hpp> // OpenSiv3D v0.6

#include "include/GraphDrawing.hpp"

class NodeVisualizer : public BasicGraphVisualizer
{
public:

	NodeVisualizer(double nodeRadius, double edgeThickness, ColorF nodeColor, ColorF edgeColor)
		: BasicGraphVisualizer(nodeRadius, edgeThickness, nodeColor, edgeColor) {}

	virtual void drawNode(const Vec2& pos, GraphEdge::IndexType) const override
	{
		const double radius = screenNodeRadius();
		const double thickness = screenFrameThickness();

		pos.asCircle(radius).draw(m_nodeColor);
		pos.asCircle(radius).drawFrame(thickness, m_edgeColor);
	}

	virtual void drawEdge(const Line&, GraphEdge::IndexType, GraphEdge::IndexType) const override {}

	double screenNodeRadius()const
	{
		return m_nodeRadius / Graphics2D::GetMaxScaling();
	}

	double screenFrameThickness()const
	{
		return m_edgeThickness / Graphics2D::GetMaxScaling();
	}
};

class EdgeVisualizer : public BasicGraphVisualizer
{
public:

	EdgeVisualizer(double edgeThickness, ColorF edgeColor)
		: BasicGraphVisualizer(0, edgeThickness, Alpha(0), edgeColor) {}

	virtual void drawNode(const Vec2&, GraphEdge::IndexType) const override {}

	virtual void drawEdge(const Line& line, GraphEdge::IndexType, GraphEdge::IndexType) const override
	{
		line.draw(screenEdgeThickness(), m_edgeColor);
	}

	double screenEdgeThickness()const
	{
		return m_edgeThickness / Graphics2D::GetMaxScaling();
	}
};

struct Theme
{
	static Color BackgroundColor() { return Color{ U"#f5f5f5" }; }

	static Color NodeColor() { return Color{ U"#eeeeee" }; }
	static Color NodeFrameColor() { return ColorF{ 0.7 }; }

	static Color EdgeColor() { return Color{ U"#a1a7b2" }; }

	static Color GradientColor(double t) { return Color{ U"#c8d5e4" }.lerp(Color{ U"#2d3743" }, t); }

	static Color AccentColor() { return Color{ U"#2692f7" }; }
};

class BFSTree
{
public:

	BFSTree() = default;

	void constructTree(const LayoutForceDirected& layout, GraphEdge::IndexType nodeIndex)
	{
		m_children.clear();

		HashSet<GraphEdge::IndexType> visited;

		m_nodeIndex = nodeIndex;

		visited.emplace(m_nodeIndex);

		for (int depthLimit = 1;; ++depthLimit)
		{
			if (!depthLimitedSearch(layout, nodeIndex, visited, depthLimit, 0))
			{
				break;
			}
		}
	}

	void drawBack(const EdgeVisualizer& edgeVis)const
	{
		drawColorEdge(edgeVis, 0);
	}

	void drawFront(const NodeVisualizer& nodeVis, const EdgeVisualizer& edgeVis, const Font& font)const
	{
		drawColorNode(nodeVis, 0);

		drawShortestPathEdge(edgeVis, 0);

		drawShortestPathNode(nodeVis, 0);

		drawShortestPathLabel(font, 0);
	}

	void turnOn(int depth = 0)
	{
		if (m_isActive && m_transitionWatch.isRunning())
		{
			return;
		}

		m_isActive = true;

		for (auto& child : m_children)
		{
			child.turnOn(depth + 1);
		}

		if (depth == 0)
		{
			m_transitionWatch.start();
		}
		else
		{
			m_transitionWatch.reset();
		}
	}

	void turnOff()
	{
		if (!m_isActive && m_transitionWatch.isRunning())
		{
			return;
		}

		m_isActive = false;

		for (auto& child : m_children)
		{
			child.turnOff();
		}

		if (m_children.empty())
		{
			m_transitionWatch.start();
		}
		else
		{
			m_transitionWatch.reset();
		}
	}

	double transitionValue()const
	{
		const double duration = 100.0;
		const double t = Math::Saturate(m_transitionWatch.msF() / duration);
		return m_isActive ? t : (1.0 - t);
	}

	void update(const HashTable<GraphEdge::IndexType, Vec2>& positions)
	{
		if (m_nodeIndex != -1)
		{
			m_position = positions.at(m_nodeIndex);
		}

		if (m_isActive)
		{
			if (1.0 <= transitionValue())
			{
				for (auto& child : m_children)
				{
					child.m_transitionWatch.start();
				}
			}
		}
		else
		{
			const bool childCompleted = m_children.all([](const BFSTree& child) { return child.transitionValue() <= 0.0; });

			if (childCompleted)
			{
				m_transitionWatch.start();
			}
		}
		
		for (auto& child : m_children)
		{
			child.update(positions);
		}
	}

	bool resetSearchIndex(Optional<GraphEdge::IndexType> searchIndex = none)
	{
		m_isOnShortestPath = searchIndex && searchIndex.value() == m_nodeIndex;

		for (auto& child : m_children)
		{
			m_isOnShortestPath = child.resetSearchIndex(searchIndex) || m_isOnShortestPath;
		}

		return m_isOnShortestPath;
	}

	GraphEdge::IndexType index()const { return m_nodeIndex; }

private:

	bool depthLimitedSearch(const LayoutForceDirected& layout, GraphEdge::IndexType nodeIndex, HashSet<GraphEdge::IndexType>& visited, int depthLimit, int depth)
	{
		bool added = false;

		if (depth < depthLimit)
		{
			for (auto p1Index : layout.activeAdjacentNodes(nodeIndex))
			{
				if (visited.find(p1Index) == visited.end())
				{
					BFSTree childTree;
					childTree.m_nodeIndex = p1Index;

					visited.emplace(p1Index);

					m_children.push_back(childTree);

					added = true;
				}
			}

			for (auto& child : m_children)
			{
				added = child.depthLimitedSearch(layout, child.m_nodeIndex, visited, depthLimit, depth + 1) || added;
			}
		}

		return added;
	}

	Color depthColor(int depth)const
	{
		const double t = Math::Saturate(1.0 * depth / 10);
		return Theme::GradientColor(t);
	}

	void drawColorNode(const NodeVisualizer& nodeVis, int depth)const
	{
		for (const auto& child : m_children)
		{
			child.drawColorNode(nodeVis, depth + 1);
		}

		m_position.asCircle(0.75 * nodeVis.screenNodeRadius() * transitionValue()).draw(depthColor(depth));
	}

	void drawColorEdge(const EdgeVisualizer& edgeVis, int depth)const
	{
		for (const auto& child : m_children)
		{
			child.drawColorEdge(edgeVis, depth + 1);
		}

		for (const auto& child : m_children)
		{
			const auto begin = m_position;
			const auto end = child.m_position;

			Line{ begin, Math::Lerp(begin, end, child.transitionValue()) }
			.draw(2 * edgeVis.screenEdgeThickness(), depthColor(depth), depthColor(depth + 1));
		}
	}

	void drawShortestPathNode(const NodeVisualizer& nodeVis, int depth)const
	{
		for (const auto& child : m_children)
		{
			child.drawShortestPathNode(nodeVis, depth + 1);
		}

		if (m_isOnShortestPath)
		{
			m_position.asCircle(nodeVis.screenNodeRadius() * transitionValue())
				.draw(Palette::White)
				.drawArc(0, 2_pi * transitionValue(), 3.0 * nodeVis.screenFrameThickness(), 0.0, Theme::AccentColor());
		}
	}

	void drawShortestPathEdge(const EdgeVisualizer& edgeVis, int depth)const
	{
		for (const auto& child : m_children)
		{
			child.drawShortestPathEdge(edgeVis, depth + 1);
		}

		for (const auto& child : m_children)
		{
			const auto begin = m_position;
			const auto end = child.m_position;

			if (child.m_isOnShortestPath)
			{
				Line{ begin , Math::Lerp(begin, end, child.transitionValue()) }
				.draw(3 * edgeVis.screenEdgeThickness(), Theme::AccentColor());
			}
		}
	}

	void drawShortestPathLabel(const Font& font, int depth)const
	{
		for (const auto& child : m_children)
		{
			child.drawShortestPathLabel(font, depth + 1);
		}

		if (m_isOnShortestPath)
		{
			const double scaling = 1.0 / Graphics2D::GetMaxScaling();

			const auto label = font(U"[", depth, U"]#", m_nodeIndex);
			const auto scale = font.fontSize() * transitionValue() * scaling;
			const Vec2 offset{ 15, -20 };

			label.draw(scale, m_position + (offset + Vec2{ 4, 4 }) * scaling, Color(U"#72706c"));
			label.draw(scale, m_position + (offset + Vec2{ 2, 2 }) * scaling, Theme::AccentColor());
			label.draw(scale, m_position + offset * scaling, Palette::White);
		}
	}

	GraphEdge::IndexType m_nodeIndex = -1;

	Vec2 m_position = Vec2{};

	Array<BFSTree> m_children;

	Stopwatch m_transitionWatch;

	bool m_isOnShortestPath = false;

	bool m_isActive = true;
};

struct GrabInfo
{
	GrabInfo() = default;

	GrabInfo(const Vec2& fromCursor, GraphEdge::IndexType index) :fromCursor(fromCursor), index(index) {}

	Vec2 fromCursor;

	GraphEdge::IndexType index;
};

void Main()
{
	Window::Resize(720, 720);

	const GraphLoader loader{ U"./clusterGraph.mtx" };

	Optional<GrabInfo> grabbingNode;

	const ForceDirectedConfig config
	{
		.autoSuspend = false,
		.initialTimeStep = 0.1,
		.updateFunction = [&](GraphEdge::IndexType nodeIndex, const Vec2& /*oldPos*/, const Vec2& newPos)
		{
			if (grabbingNode && grabbingNode.value().index == nodeIndex)
			{
				return Cursor::PosF() + grabbingNode.value().fromCursor;
			}

			return newPos;
		},
		.initialRadius = 100,
	};

	Reseed(0);

	LayoutForceDirected graph{ loader[0], config };
	graph.reset();

	Scene::SetBackground(Theme::BackgroundColor());

	NodeVisualizer nodeVisualizer{ 10.0, 1.0, Theme::NodeColor(), Theme::NodeFrameColor() };
	EdgeVisualizer edgeVisualizer{ 1.0, Theme::EdgeColor() };

	const auto drawRect = Scene::Rect().stretched(-50);

	BFSTree tree;

	Optional<GraphEdge::IndexType> startNodeIndex;
	Optional<GraphEdge::IndexType> goalNodeIndex;

	Font font(16, Typeface::Black);

	Camera2D camera{ Scene::Center(), 1.0, CameraControl::Keyboard | CameraControl::Wheel };

	Window::SetTitle(U"PathSearch -  [WASD]:カメラ操作  [マウス左]:始点を選択  [マウス右]:終点を選択  [マウス中]:ノードをつかむ");

	while (System::Update())
	{
		{
			auto t = camera.createTransformer();

			if (!MouseM.pressed())
			{
				grabbingNode = none;
				graph.setDrawArea(drawRect);
			}

			graph.update(10);

			Optional<GraphEdge::IndexType> mouseOverIndex;

			for (auto& [nodeIndex, activeNodePos] : graph.activeNodePositions())
			{
				const bool mouseOver = activeNodePos.asCircle(nodeVisualizer.m_nodeRadius / Graphics2D::GetMaxScaling()).intersects(Cursor::PosF());

				if (MouseM.down() && mouseOver)
				{
					grabbingNode = GrabInfo{ activeNodePos - Cursor::PosF(), nodeIndex };
				}

				if (!mouseOverIndex && mouseOver)
				{
					mouseOverIndex = nodeIndex;
				}
			}

			if (MouseL.down())
			{
				if (mouseOverIndex)
				{
					if (startNodeIndex != mouseOverIndex)
					{
						startNodeIndex = mouseOverIndex;
						tree.constructTree(graph, startNodeIndex.value());
						tree.turnOn();
						tree.resetSearchIndex(goalNodeIndex);
					}
				}
				else
				{
					tree.turnOff();
					startNodeIndex.reset();
					tree.resetSearchIndex();
				}
			}

			if (MouseR.down())
			{
				if (mouseOverIndex)
				{
					if (mouseOverIndex != goalNodeIndex)
					{
						goalNodeIndex = mouseOverIndex;
						tree.resetSearchIndex(goalNodeIndex);
					}
				}
				else
				{
					goalNodeIndex.reset();
					tree.resetSearchIndex();
				}
			}

			tree.update(graph.activeNodePositions());
			
			graph.draw(edgeVisualizer);

			tree.drawBack(edgeVisualizer);

			graph.draw(nodeVisualizer);

			tree.drawFront(nodeVisualizer, edgeVisualizer, font);
		}

		camera.update();
	}
}
