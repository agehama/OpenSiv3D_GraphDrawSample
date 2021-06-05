# include <Siv3D.hpp> // OpenSiv3D v0.6.0

namespace s3d
{
	struct GraphEdge
	{
		using IndexType = int32;

		IndexType index0;

		IndexType index1;

		bool operator==(const GraphEdge& e) const { return index0 == e.index0 && index1 == e.index1; }
	};

	namespace detail
	{
		template <class T>
		struct SparseEntry
		{
			SparseEntry() = default;

			SparseEntry(GraphEdge::IndexType x, GraphEdge::IndexType y, T v)
				: m_x(x)
				, m_y(y)
				, m_v(v) {}

			bool operator==(const SparseEntry& a) const { return m_x == a.m_x && m_y == a.m_y && m_v == a.m_v; }

			GraphEdge::IndexType m_x, m_y;

			T m_v;
		};

		enum class SparseFormat { CSR, CSC };

		template <class SparseEntryType>
		Array<SparseEntryType>& SortEntries(Array<SparseEntryType>& entries, SparseFormat format = SparseFormat::CSR)
		{
			switch (format)
			{
			case SparseFormat::CSR:
				std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
					return a.m_y == b.m_y ? a.m_x < b.m_x : a.m_y < b.m_y;
					});
				break;

			case SparseFormat::CSC:
				std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
					return a.m_x == b.m_x ? a.m_y < b.m_y : a.m_x < b.m_x;
					});
				break;

			default:
				break;
			}

			return entries;
		}

		template <>
		Array<GraphEdge>& SortEntries<GraphEdge>(Array<GraphEdge>& entries, SparseFormat format)
		{
			switch (format)
			{
			case SparseFormat::CSR:
				std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
					return a.index1 == b.index1 ? a.index0 < b.index0 : a.index1 < b.index1;
					});
				break;

			case SparseFormat::CSC:
				std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
					return a.index0 == b.index0 ? a.index1 < b.index1 : a.index0 < b.index0;
					});
				break;

			default:
				break;
			}

			return entries;
		}

		template <class T>
		class SparseMat
		{
		public:

			SparseMat() = default;

			explicit SparseMat(const Array<SparseEntry<T>>& sortedEntries, SparseFormat format = SparseFormat::CSR) : m_format(format)
			{
				if (sortedEntries.empty())
				{
					return;
				}

				switch (format)
				{
				case SparseFormat::CSR:
					for (const auto& entry : sortedEntries)
					{
						while (m_rowBeginIndices.size() < entry.m_y + 1)
						{
							m_rowBeginIndices.push_back(static_cast<GraphEdge::IndexType>(m_indices.size()));
						}

						m_indices.push_back(entry.m_x);
						m_values.push_back(entry.m_v);
					}
					break;

				case SparseFormat::CSC:
					for (const auto& entry : sortedEntries)
					{
						while (m_rowBeginIndices.size() < entry.m_x + 1)
						{
							m_rowBeginIndices.push_back(static_cast<GraphEdge::IndexType>(m_indices.size()));
						}

						m_indices.push_back(entry.m_y);
						m_values.push_back(entry.m_v);
					}
					break;

				default:
					break;
				}
			}

			explicit SparseMat(const Array<GraphEdge>& sortedEdges) : m_format(SparseFormat::CSR)
			{
				if (sortedEdges.empty())
				{
					return;
				}

				for (const auto& entry : sortedEdges)
				{
					while (m_rowBeginIndices.size() < entry.index1 + 1)
					{
						m_rowBeginIndices.push_back(static_cast<GraphEdge::IndexType>(m_indices.size()));
					}

					m_indices.push_back(entry.index0);
					m_values.push_back(1.f);
				}
			}

			void initFrom(const Array<SparseEntry<T>>& sortedEntries)
			{
				*this = SparseMat<T>{ sortedEntries };
			}

			void fill(GraphEdge::IndexType rowCount, GraphEdge::IndexType columnCount, T value)
			{
				if (value == 0)
				{
					m_rowBeginIndices.clear();
					m_indices.clear();
					m_values.clear();
				}
				else
				{
					const GraphEdge::IndexType elemCount = rowCount * columnCount;

					m_rowBeginIndices.resize(rowCount);
					std::iota(m_rowBeginIndices.begin(), m_rowBeginIndices.end(), 0);
					std::transform(m_rowBeginIndices.begin(), m_rowBeginIndices.end(), m_rowBeginIndices.begin(), [columnCount](GraphEdge::IndexType y) { return y * columnCount; });

					m_indices.resize(elemCount);
					std::iota(m_indices.begin(), m_indices.end(), 0);
					std::transform(m_indices.begin(), m_indices.end(), m_indices.begin(), [columnCount](GraphEdge::IndexType i) { return i % columnCount; });

					m_values.resize(elemCount);
					std::fill(m_values.begin(), m_values.end(), value);
				}
			}

			Array<SparseEntry<T>> decompressEntries()const
			{
				Array<SparseEntry<T>> entries(m_indices.size());

				switch (m_format)
				{
				case SparseFormat::CSR:
					for (GraphEdge::IndexType y = 0; y < rowCount(); ++y)
					{
						for (GraphEdge::IndexType i = rowBegin(y); i < rowEnd(y); ++i)
						{
							auto& entry = entries[i];
							entry.m_y = y;
							entry.m_x = m_indices[i];
							entry.m_v = m_values[i];
						}
					}
					break;

				case SparseFormat::CSC:
					for (GraphEdge::IndexType x = 0; x < rowCount(); ++x)
					{
						for (GraphEdge::IndexType i = rowBegin(x); i < rowEnd(x); ++i)
						{
							auto& entry = entries[i];
							entry.m_y = m_indices[i];
							entry.m_x = x;
							entry.m_v = m_values[i];
						}
					}
					break;

				default:
					break;
				}

				return entries;
			}

			SparseMat<T> multiply(SparseMat<T>& B, size_t threadCount_)
			{
				toCSR();
				B.toCSC();

				GraphEdge::IndexType* const pXA0 = &m_indices[0];
				GraphEdge::IndexType* const pXB0 = &B.m_indices[0];

				T* const pVA0 = &m_values[0];
				T* const pVB0 = &B.m_values[0];

				const auto task = [&](GraphEdge::IndexType yABegin, GraphEdge::IndexType yAEnd)
				{
					Array<SparseEntry<T>> results;

					for (GraphEdge::IndexType yA = yABegin; yA < yAEnd; ++yA)
					{
						const GraphEdge::IndexType beginIndexA = rowBegin(yA);
						const GraphEdge::IndexType endIndexA = rowEnd(yA);

						GraphEdge::IndexType* const pBeginXA = pXA0 + beginIndexA;
						GraphEdge::IndexType* const pEndXA = pXA0 + endIndexA;
						const GraphEdge::IndexType prevEndXA = *(pEndXA - 1);

						T* const pBeginVA = pVA0 + beginIndexA;

						for (GraphEdge::IndexType xB = 0; xB < B.rowCount(); ++xB)
						{
							const GraphEdge::IndexType beginIndexB = B.rowBegin(xB);
							const GraphEdge::IndexType endIndexB = B.rowEnd(xB);

							GraphEdge::IndexType* const pEndXB = pXB0 + endIndexB;
							const GraphEdge::IndexType prevEndXB = *(pEndXB - 1);

							GraphEdge::IndexType* pXA = pBeginXA;
							T* pVA = pBeginVA;

							GraphEdge::IndexType* pXB = pXB0 + beginIndexB;
							T* pVB = pVB0 + beginIndexB;

							if (pXA == pEndXA || pXB == pEndXB)
							{
								continue;
							}

							T currentValue = 0;

							for (;;)
							{
								if (*pXB < *pXA)
								{
									if (prevEndXB < *pXA || ++pXB == pEndXB)
									{
										break;
									}
									++pVB;
								}
								else if (*pXA < *pXB)
								{
									if (prevEndXA < *pXB || ++pXA == pEndXA)
									{
										break;
									}
									++pVA;
								}
								else
								{
									currentValue += *(pVA++) * *(pVB++);
									if (++pXA == pEndXA || ++pXB == pEndXB)
									{
										break;
									}
								}
							}

							if (currentValue != 0)
							{
								results.emplace_back(xB, yA, currentValue);
							}
						}
					}

					return results;
				};

				const GraphEdge::IndexType threadCount = Min(static_cast<GraphEdge::IndexType>(threadCount_), rowCount());
				const GraphEdge::IndexType rowCountPerThread = rowCount() / threadCount;

				Array<std::future<Array<SparseEntry<T>>>> futures;
				for (GraphEdge::IndexType i = 0; i < threadCount; ++i)
				{
					if (i + 1 == threadCount)
					{
						futures.push_back(std::async(std::launch::async, task, i * rowCountPerThread, rowCount()));
					}
					else
					{
						futures.push_back(std::async(std::launch::async, task, i * rowCountPerThread, (i + 1) * rowCountPerThread));
					}
				}

				Array<SparseEntry<T>> results(Arg::reserve = (rowCount() + B.rowCount()));

				for (auto& f : futures)
				{
					const auto taskResult = f.get();
					results.insert(results.end(), taskResult.begin(), taskResult.end());
				}

				return SparseMat<T>{ results, SparseFormat::CSR };
			}

			SparseMat<T> operator*(SparseMat<T>& B)
			{
				return multiply(B, 1);
			}

			SparseMat<T> operator+(const SparseMat<T>& B)const
			{
				Array<SparseEntry<T>> results;

				const GraphEdge::IndexType maxRowCount = std::max(rowCount(), B.rowCount());

				for (GraphEdge::IndexType y = 0; y < maxRowCount; ++y)
				{
					HashTable<GraphEdge::IndexType, T> currentRow;

					if (y < m_rowBeginIndices.size())
					{
						for (GraphEdge::IndexType i = rowBegin(y); i < rowEnd(y); ++i)
						{
							currentRow[m_indices[i]] = m_values[i];
						}
					}

					if (y < B.m_rowBeginIndices.size())
					{
						for (GraphEdge::IndexType i = B.rowBegin(y); i < B.rowEnd(y); ++i)
						{
							currentRow[B.m_indices[i]] += B.m_values[i];
						}
					}

					for (const auto& entry : currentRow)
					{
						results.emplace_back(entry.first, y, entry.second);
					}
				}

				return SparseMat<T>{ results, SparseFormat::CSR };
			}

			void transpose()
			{
				auto cooForm = decompressEntries();

				for (auto& entry : cooForm)
				{
					const GraphEdge::IndexType temp = entry.m_x;
					entry.m_x = entry.m_y;
					entry.m_y = temp;
				}

				*this = SparseMat{ SortEntries(cooForm, m_format), m_format };
			}

			void insert(GraphEdge::IndexType x, GraphEdge::IndexType y, T v)
			{
				auto cooForm = decompressEntries();
				cooForm.emplace_back(x, y, v);

				*this = SparseMat{ SortEntries(cooForm, m_format), m_format };
			}

			void append(const Array<SparseEntry<T>>& entries)
			{
				auto cooForm = decompressEntries();
				cooForm.insert(cooForm.end(), entries.begin(), entries.end());

				*this = SparseMat{ SortEntries(cooForm, m_format), m_format };
			}

			void toCSR()
			{
				if (m_format == SparseFormat::CSR)
				{
					return;
				}

				auto cooForm = decompressEntries();

				*this = SparseMat{ SortEntries(cooForm, SparseFormat::CSR), SparseFormat::CSR };
			}

			void toCSC()
			{
				if (m_format == SparseFormat::CSC)
				{
					return;
				}

				auto cooForm = decompressEntries();

				*this = SparseMat{ SortEntries(cooForm, SparseFormat::CSC), SparseFormat::CSC };
			}

			GraphEdge::IndexType rowCount() const noexcept { return static_cast<GraphEdge::IndexType>(m_rowBeginIndices.size()); }

			GraphEdge::IndexType rowBegin(GraphEdge::IndexType row) const { return m_rowBeginIndices[row]; }

			GraphEdge::IndexType rowEnd(GraphEdge::IndexType row) const
			{
				return m_rowBeginIndices.size() <= row + 1
					? static_cast<GraphEdge::IndexType>(m_indices.size())
					: m_rowBeginIndices[row + 1];
			}

			GraphEdge::IndexType getX(GraphEdge::IndexType i) const { return m_indices[i]; }

			T getV(GraphEdge::IndexType i) const { return m_values[i]; }

			const Array<GraphEdge::IndexType>& getRowBeginIndices() const noexcept { return m_rowBeginIndices; }

			const Array<GraphEdge::IndexType>& getXs() const noexcept { return m_indices; }

			const Array<T>& getVs() const noexcept { return m_values; }

			Array<T>& getVs() noexcept { return m_values; }

			SparseFormat getFormat() const noexcept { return m_format; }

		private:

			Array<GraphEdge::IndexType> m_rowBeginIndices;

			Array<GraphEdge::IndexType> m_indices;

			Array<T> m_values;

			SparseFormat m_format = SparseFormat::CSR;
		};

		class QuadTreeVertices
		{
		public:

			QuadTreeVertices() = default;

			explicit QuadTreeVertices(const RectF& scope)
				: m_scope(scope)
				, m_depth(0) {}

			void init(const RectF& currentScope, uint8 maxDepth, uint8 currentDepth = 0)
			{
				m_scope = currentScope;
				m_depth = currentDepth;
				m_maxDepth = maxDepth;
				m_children.clear();
				m_vertices.clear();
			}

			void add(const Vec2& vertex)
			{
				enum Dir { TL = 0, TR, BR, BL };

				const auto addToChild = [&](const Vec2& vertex)
				{
					const Vec2 center = m_scope.center();

					if (vertex.x < center.x) // Left
					{
						m_children[vertex.y < center.y ? TL : BL].add(vertex);
					}
					else // Right
					{
						m_children[vertex.y < center.y ? TR : BR].add(vertex);
					}
				};

				if (m_children) // 既に分割済みの場合
				{
					addToChild(vertex);
				}
				else if (m_vertices && (m_depth < m_maxDepth)) // 既に頂点を1つ以上持っていて最大深度より低ければ子要素に全て渡す
				{
					const Vec2 childWidth = (m_scope.size * 0.5);
					const uint8 childDepth = (m_depth + 1);

					m_children.resize(4);
					m_children[TL].init(RectF{ m_scope.pos, childWidth }, m_maxDepth, childDepth);
					m_children[TR].init(RectF{ m_scope.pos + (Vec2{ 1, 0 } *childWidth), childWidth }, m_maxDepth, childDepth);
					m_children[BR].init(RectF{ m_scope.pos + (Vec2{ 1, 1 } *childWidth), childWidth }, m_maxDepth, childDepth);
					m_children[BL].init(RectF{ m_scope.pos + (Vec2{ 0, 1 } *childWidth), childWidth }, m_maxDepth, childDepth);

					addToChild(m_vertices[0]);
					m_vertices.clear();

					addToChild(vertex);
				}
				else // 空のノードであれば頂点を1つ持つリーフノードになる
				{
					m_vertices.push_back(vertex);
				}
			}

			void update()
			{
				updateImpl();
			}

			[[nodiscard]]
			Optional<Vec2> centroid() const noexcept
			{
				if (m_numChildren == 0)
				{
					return none;
				}

				return (m_childrenSum / static_cast<double>(m_numChildren));
			}

			[[nodiscard]]
			const Vec2& size() const noexcept
			{
				return m_scope.size;
			}

			[[nodiscard]]
			bool isLeaf() const noexcept
			{
				return m_children.empty();
			}

			[[nodiscard]]
			bool isTree() const noexcept
			{
				return (not isLeaf());
			}

			[[nodiscard]]
			Vec2 repulsiveForce(const Vec2& xi, double pValue, double coeff) const
			{
				const size_t S = m_numChildren;

				if (1 <= S)
				{
					const Vec2 xS = (m_childrenSum / static_cast<double>(S));

					// Barnes-Hut opening criterion
					const double theta = 1.2;
					const double dS = Max(m_scope.w, m_scope.h);
					const bool superNode = ((dS * dS) / (xi - xS).lengthSq() <= theta * theta);

					if (superNode)
					{
						const Vec2 dist = (xS - xi);
						const double lengthSq = dist.lengthSq();

						return S * coeff * dist / Math::Pow(lengthSq, pValue * 0.5);
					}
					else if (isTree())
					{
						Vec2 f{ 0, 0 };

						for (auto& child : m_children)
						{
							f += child.repulsiveForce(xi, pValue, coeff);
						}

						return f;
					}
					else
					{
						Vec2 f{ 0, 0 };

						for (auto v : m_vertices)
						{
							if (xi == v)
							{
								continue;
							}

							const Vec2 dist = (v - xi);
							const double lengthSq = dist.lengthSq();

							f += dist / Math::Pow(lengthSq, pValue * 0.5);
						}

						return coeff * f;
					}
				}

				return { 0, 0 };
			}

		private:

			std::pair<Vec2, size_t> updateImpl()
			{
				m_childrenSum = Vec2::Zero();
				m_numChildren = 0;

				for (auto& child : m_children)
				{
					auto [v, s] = child.updateImpl();
					m_childrenSum += v;
					m_numChildren += s;
				}

				if (m_vertices)
				{
					m_childrenSum += m_vertices.sum();
					m_numChildren += m_vertices.size();
				}

				return { m_childrenSum, m_numChildren };
			}

			RectF m_scope{ 0, 0, 0, 0 };

			Vec2 m_childrenSum{ 0, 0 };

			size_t m_numChildren = 0;

			uint8 m_depth = 0;

			uint8 m_maxDepth = 0;

			Array<QuadTreeVertices> m_children;

			Array<Vec2> m_vertices;
		};

		class GraphTransform
		{
		public:

			GraphTransform() = default;

			virtual void setDrawCenter(const Vec2& drawCenter)
			{
				m_drawCenter = drawCenter;
			}

			virtual void setDrawScale(double scale)
			{
				m_drawScale = scale;
			}

			virtual void setDrawWidth(double width)
			{
				const auto rect = graphBoundingRect();
				m_drawScale = width / rect.w;
			}

			virtual void setDrawHeight(double height)
			{
				const auto rect = graphBoundingRect();
				m_drawScale = height / rect.h;
			}

			virtual void setDrawArea(const RectF& rect)
			{
				const auto graphRect = graphBoundingRect();

				m_drawScale = Min(rect.w / graphRect.w, rect.h / graphRect.h);
				m_drawCenter = rect.center();
			}

		protected:

			[[nodiscard]]
			Vec2 toDrawPos(const Vec2 graphBasePos, const Vec2& graphPos) const
			{
				return m_drawCenter + (graphPos - graphBasePos) * m_drawScale;
			}

			[[nodiscard]]
			Vec2 toGraphPos(const Vec2 graphBasePos, const Vec2& drawPos) const
			{
				return graphBasePos + (drawPos - m_drawCenter) / m_drawScale;
			}

			[[nodiscard]]
			virtual RectF graphBoundingRect() const = 0;

		private:

			// 描画位置
			Vec2 m_drawCenter{ 0, 0 };

			// 描画スケール
			double m_drawScale = 1.0;
		};

		// 連結成分ごとに分解して開始インデックスのリストを返す
		static Array<GraphEdge::IndexType> DecomposeConnectingComponents(Array<GraphEdge>& graph)
		{
			if (graph.empty())
			{
				return {};
			}

			struct UnionFind
			{
				UnionFind(GraphEdge::IndexType nodeSize)
					: parentIndex(Array<int>::IndexedGenerate(nodeSize, [](size_t i) { return static_cast<GraphEdge::IndexType>(i); }))
					, groupSize(nodeSize, 1) {}

				GraphEdge::IndexType root(GraphEdge::IndexType i)
				{
					while (parentIndex[i] != i)
					{
						i = parentIndex[i] = parentIndex[parentIndex[i]];
					}
					return i;
				}

				void unite(GraphEdge::IndexType a, GraphEdge::IndexType b)
				{
					const GraphEdge::IndexType aa = root(a);
					const GraphEdge::IndexType bb = root(b);
					if (aa == bb)
					{
						return;
					}

					if (groupSize[aa] < groupSize[bb])
					{
						groupSize[bb] += groupSize[aa];
						parentIndex[aa] = bb;
					}
					else
					{
						groupSize[aa] += groupSize[bb];
						parentIndex[bb] = aa;
					}
				}

				Array<GraphEdge::IndexType> parentIndex;

				Array<GraphEdge::IndexType> groupSize;
			};

			const auto maxIt = std::max_element(graph.begin(), graph.end(), [](const auto& a, const auto& b) { return Max(a.index0, a.index1) < Max(b.index0, b.index1); });

			const GraphEdge::IndexType nodeSize = Max(maxIt->index0, maxIt->index1) + 1;

			UnionFind unionFind(nodeSize);

			for (const auto& edge : graph)
			{
				unionFind.unite(edge.index0, edge.index1);
			}

			std::sort(graph.begin(), graph.end(), [&](const auto& a, const auto& b) { return unionFind.root(a.index0) < unionFind.root(b.index0); });

			HashSet<GraphEdge::IndexType> added;
			Array<GraphEdge::IndexType> beginIndices;

			size_t i = 0;
			for (; i < graph.size(); ++i)
			{
				const GraphEdge::IndexType rootIndex = unionFind.root(graph[i].index0);

				if (added.find(rootIndex) == added.end())
				{
					beginIndices.push_back(static_cast<GraphEdge::IndexType>(i));
					added.insert(rootIndex);
				}
			}

			beginIndices.push_back(static_cast<GraphEdge::IndexType>(i));

			return beginIndices;
		}

		// グラフから参照されていないインデックスを詰めてエッジリストとノード数を返す
		static std::pair<Array<GraphEdge>, size_t> ShrinkVertices(typename Array<GraphEdge>::const_iterator begin, typename Array<GraphEdge>::const_iterator end)
		{
			if (begin == end)
			{
				return {};
			}

			GraphEdge::IndexType maxVertexIndex = 0;

			HashSet<GraphEdge::IndexType> referencedVertices;

			for (auto it = begin; it != end; ++it)
			{
				const auto& edge = *it;

				const GraphEdge::IndexType index0 = static_cast<GraphEdge::IndexType>(edge.index0);
				const GraphEdge::IndexType index1 = static_cast<GraphEdge::IndexType>(edge.index1);

				referencedVertices.insert(index0);
				referencedVertices.insert(index1);

				maxVertexIndex = Max({ index0, index1, maxVertexIndex });
			}

			GraphEdge::IndexType skippedVertexCount = 0;

			HashTable<GraphEdge::IndexType, GraphEdge::IndexType> indexReplaceMap;

			for (GraphEdge::IndexType i = 0; i <= maxVertexIndex; ++i)
			{
				if (referencedVertices.find(i) == referencedVertices.end())
				{
					++skippedVertexCount;
				}
				else
				{
					indexReplaceMap[i] = i - skippedVertexCount;
				}
			}

			Array<GraphEdge> newGraph(Arg::reserve = std::distance(begin, end));

			for (auto it = begin; it != end; ++it)
			{
				auto newRowIt = indexReplaceMap.find(it->index1);
				auto newColIt = indexReplaceMap.find(it->index0);

				newGraph.emplace_back(newRowIt->second, newColIt->second);
			}

			return{ newGraph, indexReplaceMap.size() };
		}
	}

	class ConnectedGraph
	{
	public:

		ConnectedGraph() = default;

		ConnectedGraph(const Array<GraphEdge>& edgeList)
			: ConnectedGraph(edgeList.begin(), edgeList.end())
		{}

		template <class ForwardIt>
		ConnectedGraph(ForwardIt first, ForwardIt last)
		{
			Array<GraphEdge> edges(Arg::reserve = (std::distance(first, last) * 2));

			for (ForwardIt it = first; it != last; ++it)
			{
				edges.emplace_back(it->index0, it->index1);
				edges.emplace_back(it->index1, it->index0);
			}

			detail::SortEntries(edges);

			edges.unique_consecutive();

			Array<GraphEdge::IndexType> componentIndices = detail::DecomposeConnectingComponents(edges);

			Array<GraphEdge::IndexType> componentEdgeCounts;

			std::adjacent_difference(componentIndices.begin(), componentIndices.end(), std::back_inserter(componentEdgeCounts));

			{
				const GraphEdge::IndexType beginIndex = componentIndices[0];
				const GraphEdge::IndexType componentEdgeCount = componentEdgeCounts[1];

				std::tie(m_edges, m_nodeCount) = detail::ShrinkVertices(edges.cbegin() + beginIndex, edges.cbegin() + beginIndex + componentEdgeCount);

				detail::SortEntries(m_edges);
			}
		}

		const Array<GraphEdge>& edges() const
		{
			return m_edges;
		}

		size_t nodeCount() const
		{
			return m_nodeCount;
		}

	private:

		ConnectedGraph(const Array<GraphEdge>& edgeList, size_t nodeCount)
			: m_edges(edgeList)
			, m_nodeCount(nodeCount) {}

		friend class GraphLoader;

		Array<GraphEdge> m_edges;

		size_t m_nodeCount;
	};

	struct IGraphVisualizer
	{
		virtual ~IGraphVisualizer() = default;

		virtual void drawNode(const Vec2& pos, GraphEdge::IndexType nodeIndex) const = 0;

		virtual void drawEdge(const Line& line, GraphEdge::IndexType nodeIndex0, GraphEdge::IndexType nodeIndex1) const = 0;
	};

	class BasicGraphVisualizer : public IGraphVisualizer
	{
	public:

		explicit BasicGraphVisualizer(double nodeRadius = 10.0, double edgeThickness = 1.0, ColorF nodeColor = Palette::White, ColorF edgeColor = ColorF(0.8))
			: m_nodeRadius(nodeRadius)
			, m_edgeThickness(edgeThickness)
			, m_nodeColor(nodeColor)
			, m_edgeColor(edgeColor) {}

		virtual ~BasicGraphVisualizer() = default;

		virtual void drawNode(const Vec2& pos, GraphEdge::IndexType) const override
		{
			pos.asCircle(m_nodeRadius).draw(m_nodeColor);
		}

		virtual void drawEdge(const Line& line, GraphEdge::IndexType, GraphEdge::IndexType) const override
		{
			line.draw(m_edgeThickness, m_edgeColor);
		}

		double m_nodeRadius;

		double m_edgeThickness;

		ColorF m_nodeColor;

		ColorF m_edgeColor;
	};

	class GraphLoader
	{
	public:

		GraphLoader() = default;

		explicit GraphLoader(const FilePath& path)
		{
			if (FileSystem::Extension(path) == U"mtx")
			{
				loadMMCoordinateFormat(path);
			}
			else
			{
				loadEdgeList(path);
			}
		}

		explicit GraphLoader(const Array<GraphEdge>& edgeList)
		{
			loadEdgeList(edgeList.begin(), edgeList.end());
		}

		template <class ForwardIt>
		GraphLoader(ForwardIt first, ForwardIt last)
		{
			loadEdgeList(first, last);
		}

		// .mtx (Matrix Market Exchange Formats) からロード
		void loadMMCoordinateFormat(const FilePath& path)
		{
			TextReader reader{ path };

			if (not reader)
			{
				throw Error{ U"LoadMMCoordinateFormat(): Failed to open `{0}`"_fmt(path) };
			}

			size_t edgeCount = 0;
			{
				String line;

				while (reader.readLine(line))
				{
					if (line.starts_with(U"%")) // comment
					{
						continue;
					}

					const auto nums = line.split(U' ').map([](const String& str) { return ParseInt<int32>(str); });

					if (nums.size() == 3)
					{
						edgeCount = nums[2];
						break;
					}
				}
			}

			Array<GraphEdge> edges(Arg::reserve = edgeCount);
			{
				String line;

				while (reader.readLine(line))
				{
					if (line.starts_with(U"%")) // comment
					{
						continue;
					}

					const auto numStrs = line.split(U' ');
					if (2 <= numStrs.size())
					{
						const GraphEdge::IndexType num0 = ParseInt<GraphEdge::IndexType>(numStrs[0]);
						const GraphEdge::IndexType num1 = ParseInt<GraphEdge::IndexType>(numStrs[1]);

						if (num0 != num1)
						{
							edges.emplace_back(num0, num1);
							edges.emplace_back(num1, num0);
						}
					}
				}
			}

			loadDecomposed(edges);
		}

		// エッジのペアのリストをロード
		void loadEdgeList(const FilePath& path)
		{
			TextReader reader{ path };

			if (not reader)
			{
				throw Error{ U"loadEdgeList(): Failed to open `{0}`"_fmt(path) };
			}

			Array<GraphEdge> edges;

			String line;

			while (reader.readLine(line))
			{
				const auto elements = line.split(U' ').map([](const String& str) { return Parse<GraphEdge::IndexType>(str); });

				if (2 <= elements.size())
				{
					const GraphEdge::IndexType indexA = elements[0];
					const GraphEdge::IndexType indexB = elements[1];

					if (indexA != indexB)
					{
						edges.emplace_back(indexA, indexB);
						edges.emplace_back(indexB, indexA);
					}
				}
			}

			loadDecomposed(edges);
		}

		template <class ForwardIt>
		void loadEdgeList(ForwardIt first, ForwardIt last)
		{
			Array<GraphEdge> edges(Arg::reserve = (std::distance(first, last) * 2));

			for (ForwardIt it = first; it != last; ++it)
			{
				edges.emplace_back(it->index0, it->index1);
				edges.emplace_back(it->index1, it->index0);
			}

			loadDecomposed(edges);
		}

		size_t connectedComponentCount() const noexcept
		{
			return m_connectedComponents.size();
		}

		size_t nodeCount(size_t connectedComponentIndex) const
		{
			return m_connectedComponents[connectedComponentIndex].nodeCount();
		}

		const ConnectedGraph& connectedComponent(size_t i) const
		{
			return m_connectedComponents[i];
		}

		ConnectedGraph& connectedComponent(size_t i)
		{
			return m_connectedComponents[i];
		}

		const ConnectedGraph& operator[](size_t i)const
		{
			return m_connectedComponents[i];
		}

		size_t size() const noexcept
		{
			return m_connectedComponents.size();
		}

		Array<ConnectedGraph>::iterator begin() noexcept
		{
			return m_connectedComponents.begin();
		}

		Array<ConnectedGraph>::const_iterator begin() const noexcept
		{
			return m_connectedComponents.begin();
		}

		Array<ConnectedGraph>::iterator end() noexcept
		{
			return m_connectedComponents.end();
		}

		Array<ConnectedGraph>::const_iterator end() const noexcept
		{
			return m_connectedComponents.end();
		}

	private:

		void loadDecomposed(Array<GraphEdge>& edges)
		{
			Array<GraphEdge::IndexType> componentIndices = detail::DecomposeConnectingComponents(edges);

			Array<GraphEdge::IndexType> componentEdgeCounts;

			std::adjacent_difference(componentIndices.begin(), componentIndices.end(), std::back_inserter(componentEdgeCounts));

			for (size_t i = 0; i + 1 < componentIndices.size(); ++i)
			{
				const GraphEdge::IndexType beginIndex = componentIndices[i];
				const GraphEdge::IndexType componentEdgeCount = componentEdgeCounts[i + 1];

				auto result = detail::ShrinkVertices(edges.cbegin() + beginIndex, edges.cbegin() + beginIndex + componentEdgeCount);

				detail::SortEntries(result.first);

				m_connectedComponents.push_back(ConnectedGraph(std::move(result.first), result.second));
			}
		}

		Array<ConnectedGraph> m_connectedComponents;
	};

	struct ForceDirectedConfig
	{
		// ノード同士の斥力がどれだけ遠くまで働くかを決める値
		//  小さいと遠距離間のノードも影響し合う [2.0, 3.0] くらい
		double repulsiveExponent = 2.5;

		// 収束判定に用いる係数
		//  小さいほうが厳密になるがその分時間がかかる
		//  0.01 <= tol <= 1.0 くらい
		double tol = 0.1;

		// タイムステップを内部で調整する時のスケール
		//  小さいほど大雑把に計算を行うため速くなるが小さすぎると振動して収束しなくなる
		//  ずっと振動してなかなか収束しない場合は1に近づけると改善するかもしれない
		//  0.9 <= stepScale <= 0.99 くらい
		double stepScale = 0.9;

		// 計算に用いるスレッド数
		//  none だと内部で自動的に決める
		//  ノード数1000未満くらいの規模なら1スレッドの方が速い
		//  大きいときはスレッド数を増やすほど速くなる
		Optional<size_t> threadCount = none;

		// 自動停止フラグ
		//  true だと収束条件を満たしたらシミュレーションを終了する
		//  完了後も頂点の移動を行う場合は false に設定する
		bool autoSuspend = true;

		// 初期タイムステップ幅
		//  最適なレイアウトになる前に収束する場合はこの値を増やす
		//  動きが大きすぎてレイアウトが崩れる場合はこの値を減らす
		//  0.01 <= initialTimeStep <= 10.0 くらい
		//  ノードに外力を加えるようなインタラクティブなアプリでは小さめ [0.01, 0.1] くらいが推奨
		double initialTimeStep = 1.0;

		// ノード座標の更新関数:
		// 	設定するとノード位置の更新を上書きできる
		//  第3引数の pos をそのまま返すと何もしないのと同じ
		std::function<Vec2(GraphEdge::IndexType /*nodeIndex*/, const Vec2& /*oldPos*/, const Vec2& /*pos*/)> updateFunction;

		// ↓は基本変える必要ない
		double initialRadius = 2500.0;
		double C = 0.2;
		double K = 1.0;
		uint8 maxQuadTreeDepth = 8;
	};

	class LayoutForceDirected : public detail::GraphTransform
	{
	public:

		LayoutForceDirected() = default;

		template <class URBG>
		LayoutForceDirected(const ConnectedGraph& connectedGraph, URBG&& urbg, const ForceDirectedConfig& config)
		{
			init(connectedGraph, urbg, config);
		}

		LayoutForceDirected(const ConnectedGraph& connectedGraph, const ForceDirectedConfig& config)
		{
			init(connectedGraph, GetDefaultRNG(), config);
		}

		template <class URBG>
		void init(const ConnectedGraph& connectedGraph, URBG&& urbg, const ForceDirectedConfig& config)
		{
			m_config = config;

			Array<detail::SparseEntry<float>> points(Arg::reserve = (connectedGraph.nodeCount() * 2));

			for (auto i : step(static_cast<GraphEdge::IndexType>(connectedGraph.nodeCount())))
			{
				const auto v = static_cast<Vector2D<float>>(RandomVec2(Circle{ Vec2::Zero(), m_config.initialRadius }, std::forward<URBG>(urbg)));
				points.emplace_back(0, i, v.x);
				points.emplace_back(1, i, v.y);
			}

			m_positions = detail::SparseMat<float>{ points };
			m_oldPositions = m_positions;

			m_adjacencyMatrix = detail::SparseMat<float>{ connectedGraph.edges() };

			m_originalNodeIndices.resize(connectedGraph.nodeCount());
			std::iota(m_originalNodeIndices.begin(), m_originalNodeIndices.end(), 0);

			makeCoarseGraphSeries(std::forward<URBG>(urbg));
		}

		template <class URBG>
		bool update(int32 maxUpdateMillisec, URBG&& urbg)
		{
			if (m_converged && (1 <= m_graphDepth || 0 == m_graphDepth && m_config.autoSuspend))
			{
				return true;
			}

			Stopwatch watch{ StartImmediately::Yes };

			if (m_coarserGraph && (not m_coarserGraph->m_converged))
			{
				m_coarserGraph->update(maxUpdateMillisec);

				if (m_coarserGraph->m_converged)
				{
					copyFromCoarserGraph(std::forward<URBG>(urbg));
					prepareSimulation();
				}
			}

			if ((not m_coarserGraph) || m_coarserGraph->m_converged)
			{
				while (maxUpdateMillisec <= 0 || watch.ms() < maxUpdateMillisec)
				{
					simulateOneStep();

					if (m_config.autoSuspend && m_converged || maxUpdateMillisec == 0)
					{
						break;
					}
				}
			}

			m_elapsedSec += watch.sF();

			return m_converged;
		}

		// maxUpdateMillisec を上回るまで計算を続ける
		// 0を指定すると1回だけ計算を行う
		// -1を指定するとレイアウトが完了するまで計算を続ける
		bool update(int32 maxUpdateMillisec = 10)
		{
			return update(maxUpdateMillisec, GetDefaultRNG());
		}

		[[nodiscard]]
		GraphEdge::IndexType nodeCount() const noexcept
		{
			return m_positions.rowCount();
		}

		[[nodiscard]]
		GraphEdge::IndexType edgeCount() const noexcept
		{
			return static_cast<GraphEdge::IndexType>(m_adjacencyMatrix.getXs().size());
		}

		// [0.0, 1.0] の範囲でレイアウトの進み具合を返す
		[[nodiscard]]
		double progressRate() const noexcept
		{
			if (isActiveDepth())
			{
				// 最後は1.0を少し超える
				return Min(1.0, m_currentProgressRate);
			}

			return m_coarserGraph->progressRate();
		}

		[[nodiscard]]
		bool isCompleted() const noexcept
		{
			return m_converged;
		}

		// 計算に失敗したらtrueを返す（レイアウトの品質とは無関係）
		// 複雑なグラフだと乱数やconfigによっては振動を続けて収束しなくなることがあるため、一定期間進まなかったら失敗とみなしている
		[[nodiscard]]
		bool isFailed() const noexcept
		{
			if (m_failed)
			{
				return true;
			}
			else if (m_converged)
			{
				return false;
			}

			if (m_coarserGraph)
			{
				return m_coarserGraph->isFailed();
			}

			return false;
		}

		[[nodiscard]]
		double simulationTime() const noexcept
		{
			return m_elapsedSec;
		}

		// レイアウトの計算は少ないノード数から始まり実行が進むごとに増えていく
		// 現在処理中のグラフの情報を active~~() で取得する

		[[nodiscard]]
		Array<Vec2> activeNodePositions() const
		{
			if (isActiveDepth())
			{
				const Vec2 graphCenter = graphBoundingRect().center();

				Array<Vec2> points(nodeCount());

				for (auto [nodeIndex, position] : IndexedRef(points))
				{
					const GraphEdge::IndexType internalIndex = m_positions.rowBegin(static_cast<GraphEdge::IndexType>(nodeIndex));
					const Vec2 pos{ m_positions.getV(internalIndex), m_positions.getV(internalIndex + 1) };

					position = toDrawPos(graphCenter, pos);
				}

				return points;
			}

			return m_coarserGraph->activeNodePositions();
		}

		[[nodiscard]]
		Array<std::pair<GraphEdge::IndexType, GraphEdge::IndexType>> activeEdges() const
		{
			if (isActiveDepth())
			{
				Array<std::pair<GraphEdge::IndexType, GraphEdge::IndexType>> edges(Arg::reserve = edgeCount());

				for (GraphEdge::IndexType p0Index = 0; p0Index < nodeCount(); ++p0Index)
				{
					for (GraphEdge::IndexType i = m_adjacencyMatrix.rowBegin(p0Index); i < m_adjacencyMatrix.rowEnd(p0Index); ++i)
					{
						const GraphEdge::IndexType p1Index = m_adjacencyMatrix.getX(i);
						if (p0Index <= p1Index)
						{
							continue;
						}

						edges.emplace_back(p0Index, p1Index);
					}
				}

				return edges;
			}

			return m_coarserGraph->activeEdges();
		}

		[[nodiscard]]
		Array<GraphEdge::IndexType> activeAdjacentNodes(GraphEdge::IndexType activeNodeIndex) const
		{
			if (isActiveDepth())
			{
				return adjacentNodes(activeNodeIndex);
			}

			return m_coarserGraph->activeAdjacentNodes(activeNodeIndex);
		}

		// 現在処理中のノードインデックスを元のグラフの対応するインデックスに変換する
		[[nodiscard]]
		GraphEdge::IndexType originalNodeIndex(GraphEdge::IndexType activeNodeIndex) const
		{
			if (isActiveDepth())
			{
				return m_originalNodeIndices[activeNodeIndex];
			}

			return m_coarserGraph->originalNodeIndex(activeNodeIndex);
		}


		void draw(const IGraphVisualizer& visualizer) const
		{
			if (isActiveDepth())
			{
				const Vec2 graphCenter = graphBoundingRect().center();

				for (GraphEdge::IndexType p0Index = 0; p0Index < nodeCount(); ++p0Index)
				{
					for (GraphEdge::IndexType i = m_adjacencyMatrix.rowBegin(p0Index); i < m_adjacencyMatrix.rowEnd(p0Index); ++i)
					{
						const GraphEdge::IndexType p1Index = m_adjacencyMatrix.getX(i);
						if (p0Index <= p1Index)
						{
							continue;
						}

						const GraphEdge::IndexType p0InternalIndex = m_positions.rowBegin(p0Index);
						const GraphEdge::IndexType p1InternalIndex = m_positions.rowBegin(p1Index);

						const Vec2 p0{ m_positions.getV(p0InternalIndex), m_positions.getV(p0InternalIndex + 1) };
						const Vec2 p1{ m_positions.getV(p1InternalIndex), m_positions.getV(p1InternalIndex + 1) };

						visualizer.drawEdge(Line{ toDrawPos(graphCenter, p0), toDrawPos(graphCenter, p1) }, m_originalNodeIndices[p0Index], m_originalNodeIndices[p1Index]);
					}
				}

				for (GraphEdge::IndexType nodeIndex = 0; nodeIndex < nodeCount(); ++nodeIndex)
				{
					const GraphEdge::IndexType internalIndex = m_positions.rowBegin(nodeIndex);
					const Vec2 pos{ m_positions.getV(internalIndex), m_positions.getV(internalIndex + 1) };

					visualizer.drawNode(toDrawPos(graphCenter, pos), m_originalNodeIndices[nodeIndex]);
				}
			}

			if (m_coarserGraph)
			{
				m_coarserGraph->draw(visualizer);
			}
		}

		[[nodiscard]]
		RectF boundingRect() const
		{
			const auto graphRect = graphBoundingRect();
			const Vec2 tl = toDrawPos(graphRect.center(), graphRect.tl());
			const Vec2 br = toDrawPos(graphRect.center(), graphRect.br());
			return RectF{ tl, br - tl };
		}

		[[nodiscard]]
		Vec2 centroid() const
		{
			const auto graphRect = graphBoundingRect();
			const auto pos = graphCentroid();
			return toDrawPos(graphRect.center(), pos);
		}

		template <class URBG>
		void reset(URBG&& urbg)
		{
			auto& vs = m_positions.getVs();
			for (auto i : step(nodeCount()))
			{
				const auto v = static_cast<Vector2D<float>>(RandomVec2(Circle{ Vec2::Zero(), m_config.initialRadius }, std::forward<URBG>(urbg)));
				const GraphEdge::IndexType internalIndex = m_positions.rowBegin(i);
				vs[internalIndex + 0] = v.x;
				vs[internalIndex + 1] = v.y;
			}

			resetImpl();
		}

		void reset()
		{
			reset(GetDefaultRNG());
		}

		void setDrawCenter(const Vec2& drawCenter) override
		{
			GraphTransform::setDrawCenter(drawCenter);

			if (m_coarserGraph)
			{
				m_coarserGraph->setDrawCenter(drawCenter);
			}
		}

		void setDrawScale(double scale) override
		{
			GraphTransform::setDrawScale(scale);

			if (m_coarserGraph)
			{
				m_coarserGraph->setDrawScale(scale);
			}
		}

		void setDrawWidth(double width) override
		{
			GraphTransform::setDrawWidth(width);

			if (m_coarserGraph)
			{
				m_coarserGraph->setDrawWidth(width);
			}
		}

		void setDrawHeight(double height) override
		{
			GraphTransform::setDrawHeight(height);

			if (m_coarserGraph)
			{
				m_coarserGraph->setDrawHeight(height);
			}
		}

		void setDrawArea(const RectF& rect) override
		{
			GraphTransform::setDrawArea(rect);

			if (m_coarserGraph)
			{
				m_coarserGraph->setDrawArea(rect);
			}
		}

	private:

		[[nodiscard]]
		RectF graphBoundingRect() const override
		{
			if (nodeCount() == 0)
			{
				return Rect{};
			}

			if (isActiveDepth())
			{
				float minX = Largest<float>;
				float minY = Largest<float>;
				float maxX = -Largest<float>;
				float maxY = -Largest<float>;

				for (GraphEdge::IndexType nodeIndex = 0; nodeIndex < nodeCount(); ++nodeIndex)
				{
					const GraphEdge::IndexType internalIndex = m_positions.rowBegin(nodeIndex);
					const float x = m_positions.getV(internalIndex + 0);
					const float y = m_positions.getV(internalIndex + 1);

					minX = Min(minX, x);
					minY = Min(minY, y);
					maxX = Max(maxX, x);
					maxY = Max(maxY, y);
				}

				const Vec2 minPos{ minX, minY };
				const Vec2 maxPos{ maxX, maxY };
				return RectF{ minPos, maxPos - minPos };
			}

			if (m_coarserGraph)
			{
				return m_coarserGraph->graphBoundingRect();
			}

			return Rect{};
		}

		[[nodiscard]]
		Vec2 graphCentroid() const
		{
			if (nodeCount() == 0)
			{
				return Vec2{};
			}

			if (isActiveDepth())
			{
				Vec2 sum{ 0, 0 };

				for (GraphEdge::IndexType nodeIndex = 0; nodeIndex < nodeCount(); ++nodeIndex)
				{
					const GraphEdge::IndexType internalIndex = m_positions.rowBegin(nodeIndex);
					const float x = m_positions.getV(internalIndex + 0);
					const float y = m_positions.getV(internalIndex + 1);
					sum += Vec2{ x, y };
				}

				return sum / nodeCount();
			}

			if (m_coarserGraph)
			{
				return m_coarserGraph->graphCentroid();
			}

			return Vec2{};
		}

		[[nodiscard]]
		bool isActiveDepth() const
		{
			/*
			アクティブなグラフは次のいずれか
			- 自分が収束していなくて子供がいない（初期状態）
			- 自分が収束していなくて子供が収束している（計算中）
			- 自分が収束していて一番上の親（完了後）
			*/
			return (not m_converged) && ((not m_coarserGraph) || m_coarserGraph && m_coarserGraph->m_converged)
				|| m_converged && m_graphDepth == 0;
		}

		[[nodiscard]]
		size_t useThreadCount() const noexcept
		{
			if (m_config.threadCount)
			{
				return m_config.threadCount.value();
			}

			// 1000ノード以下ならシングルスレッドの方が速い
			return nodeCount() < 1000 ? 1 : std::thread::hardware_concurrency();
		}

		[[nodiscard]]
		Array<GraphEdge::IndexType> adjacentNodes(GraphEdge::IndexType nodeIndex) const
		{
			Array<GraphEdge::IndexType> indices;

			for (GraphEdge::IndexType i = m_adjacencyMatrix.rowBegin(nodeIndex); i < m_adjacencyMatrix.rowEnd(nodeIndex); ++i)
			{
				indices.push_back(m_adjacencyMatrix.getX(i));
			}

			return indices;
		}

		void resetImpl()
		{
			m_elapsedSec = 0.0;
			m_converged = false;
			m_failed = false;

			// シードごとに同じ結果になるようにinit()で初期化したときと乱数の呼び出し回数を合わせる
			for (size_t i = 1; i < m_initRandomCount; ++i)
			{
				const double volatile t = Random();
			}

			if (m_coarserGraph)
			{
				auto P_T = m_prolongationMatrix;
				P_T.transpose();
				m_coarserGraph->m_positions = P_T * m_positions;
				m_positions.toCSR();

				auto& vs = m_coarserGraph->m_positions.getVs();
				const auto& scaleVs = m_coarserGraph->m_scaleVector.getVs();
				for (auto i : step(m_coarserGraph->nodeCount()))
				{
					const GraphEdge::IndexType internalIndex = m_coarserGraph->m_positions.rowBegin(i);
					vs[internalIndex + 0] /= scaleVs[i];
					vs[internalIndex + 1] /= scaleVs[i];
				}

				m_coarserGraph->resetImpl();
			}
		}

		/*
		ここからinit()で呼ばれるデータ構築用の関数
		*/

		template <class URBG>
		void makeCoarseGraphSeries(URBG&& urbg)
		{
			coarsenByEC(std::forward<URBG>(urbg));

			const double threshold = 0.75;
			if (1.0 * m_coarserGraph->nodeCount() / nodeCount() > threshold)
			{
				m_prolongationMatrix.fill(0, 0, 0);

				coarsenByMIVS();
			}

			calculateAdjacentMat();

			const double rate = 1.0 * m_coarserGraph->nodeCount() / nodeCount();

			m_coarserGraph->m_config = m_config;

			const double gamma = Math::Sqrt(7.0 / 4.0);
			m_coarserGraph->m_config.K = m_config.K * gamma;

			const double cutoff = 0.5;
			if (cutoff < rate)
			{
				const double gamma_ = Math::Lerp(gamma, 1.0, Math::InvLerp(cutoff, 1.0, rate));
				m_coarserGraph->m_config.K = m_config.K * gamma_;
			}

			m_coarserGraph->m_graphDepth = m_graphDepth + 1;

			const GraphEdge::IndexType minNodeCount = 3;
			if (m_coarserGraph->nodeCount() < minNodeCount || m_coarserGraph->nodeCount() == nodeCount())
			{
				m_coarserGraph.reset();
			}
			else if (minNodeCount < m_coarserGraph->nodeCount())
			{
				m_coarserGraph->makeCoarseGraphSeries(std::forward<URBG>(urbg));
			}
		}

		template <class URBG>
		void coarsenByEC(URBG&& urbg)
		{
			m_coarserGraph = std::make_unique<LayoutForceDirected>();

			Array<std::pair<GraphEdge::IndexType, GraphEdge::IndexType>> maximalMatching;
			HashSet<GraphEdge::IndexType> matchedVertices;
			HashSet<GraphEdge::IndexType> unmatchedVertices;
			for (GraphEdge::IndexType y = 0; y < nodeCount(); ++y)
			{
				unmatchedVertices.emplace(y);
			}

			Array<GraphEdge::IndexType> visitingVertices(nodeCount());
			std::iota(visitingVertices.begin(), visitingVertices.end(), 0);
			visitingVertices.shuffle(std::forward<URBG>(urbg));
			m_initRandomCount = visitingVertices.size();

			struct RowVal
			{
				RowVal(GraphEdge::IndexType x, float val) :x(x), val(val) {}
				GraphEdge::IndexType x;
				float val;
			};

			for (GraphEdge::IndexType nodeIndex : visitingVertices)
			{
				if (matchedVertices.end() != matchedVertices.find(nodeIndex))
				{
					continue;
				}

				Array<RowVal> currentRow;

				for (GraphEdge::IndexType i = m_adjacencyMatrix.rowBegin(nodeIndex); i < m_adjacencyMatrix.rowEnd(nodeIndex); ++i)
				{
					currentRow.emplace_back(m_adjacencyMatrix.getX(i), m_adjacencyMatrix.getV(i));
				}

				std::sort(currentRow.begin(), currentRow.end(), [](const auto& a, const auto& b) { return a.val > b.val; });

				for (const auto& elem : currentRow)
				{
					if (matchedVertices.end() == matchedVertices.find(elem.x))
					{
						matchedVertices.emplace(elem.x);
						matchedVertices.emplace(nodeIndex);

						unmatchedVertices.erase(elem.x);
						unmatchedVertices.erase(nodeIndex);

						maximalMatching.emplace_back(elem.x, nodeIndex);

						break;
					}
				}
			}

			for (auto fineIndex : unmatchedVertices)
			{
				maximalMatching.emplace_back(fineIndex, fineIndex);
			}

			Array<detail::SparseEntry<float>> pEntries(Arg::reserve = (m_adjacencyMatrix.rowCount() * maximalMatching.size()));

			m_coarserGraph->m_originalNodeIndices.resize(maximalMatching.size());

			for (const auto& [coarserIndex, edge] : Indexed(maximalMatching))
			{
				pEntries.emplace_back(static_cast<GraphEdge::IndexType>(coarserIndex), edge.first, 1.0f);
				if (edge.first != edge.second)
				{
					pEntries.emplace_back(static_cast<GraphEdge::IndexType>(coarserIndex), edge.second, 1.0f);
				}

				m_coarserGraph->m_originalNodeIndices[coarserIndex] = m_originalNodeIndices[edge.first];
			}

			m_prolongationMatrix = detail::SparseMat<float>{ SortEntries(pEntries) };
		}

		void coarsenByMIVS()
		{
			m_coarserGraph = std::make_unique<LayoutForceDirected>();

			HashTable<GraphEdge::IndexType, GraphEdge::IndexType> maxIndependentVertices;

			// fine vertex id -> depend degreees
			Array<int32> dependDegrees(nodeCount(), 0);

			// fine vertex id -> [child fine vertex id]
			Array<HashSet<GraphEdge::IndexType>> children(nodeCount(), HashSet<GraphEdge::IndexType>{});

			for (GraphEdge::IndexType y = 0; y < nodeCount(); ++y)
			{
				if (dependDegrees[y] == 0)
				{
					const GraphEdge::IndexType coarseIndex = static_cast<GraphEdge::IndexType>(maxIndependentVertices.size());
					maxIndependentVertices.emplace(y, coarseIndex);
					dependDegrees[y] = 1;

					for (GraphEdge::IndexType i = m_adjacencyMatrix.rowBegin(y); i < m_adjacencyMatrix.rowEnd(y); ++i)
					{
						const GraphEdge::IndexType x = m_adjacencyMatrix.getX(i);

						if (maxIndependentVertices.end() == maxIndependentVertices.find(x))
						{
							++dependDegrees[x];
							children[y].emplace(x);
						}
					}
				}
			}

			Array<detail::SparseEntry<float>> pEntries(Arg::reserve = (m_adjacencyMatrix.rowCount() * maxIndependentVertices.size()));

			m_coarserGraph->m_originalNodeIndices.resize(maxIndependentVertices.size());

			for (const auto [finerIndex, coarserIndex] : maxIndependentVertices)
			{
				pEntries.emplace_back(coarserIndex, finerIndex, 1.0f);

				for (auto childIndex : children[finerIndex])
				{
					const int32 childDegrees = dependDegrees[childIndex];
					pEntries.emplace_back(coarserIndex, childIndex, 1.0f / childDegrees);
				}

				m_coarserGraph->m_originalNodeIndices[coarserIndex] = m_originalNodeIndices[finerIndex];
			}

			m_prolongationMatrix = detail::SparseMat<float>{ SortEntries(pEntries) };
		}

		void calculateAdjacentMat()
		{
			auto transposedProlongationMatrix = m_prolongationMatrix;
			transposedProlongationMatrix.transpose();

			m_coarserGraph->m_adjacencyMatrix = transposedProlongationMatrix.multiply(m_adjacencyMatrix, useThreadCount()).multiply(m_prolongationMatrix, useThreadCount());

			detail::SparseMat<float> scaleVector;
			scaleVector.fill(nodeCount(), 1, 1);

			m_coarserGraph->m_scaleVector = transposedProlongationMatrix * scaleVector;

			m_coarserGraph->m_positions = transposedProlongationMatrix * m_positions;
			m_positions.toCSR();

			auto& vs = m_coarserGraph->m_positions.getVs();
			const auto& scaleVs = m_coarserGraph->m_scaleVector.getVs();
			for (auto i : step(m_coarserGraph->nodeCount()))
			{
				const GraphEdge::IndexType internalIndex = m_coarserGraph->m_positions.rowBegin(i);
				vs[internalIndex + 0] /= scaleVs[i];
				vs[internalIndex + 1] /= scaleVs[i];
			}

			m_coarserGraph->m_oldPositions = m_coarserGraph->m_positions;

			m_adjacencyMatrix.toCSR();

			// 対角成分を0で埋める
			auto& coarserAdjacencyMatrix = m_coarserGraph->m_adjacencyMatrix;
			for (GraphEdge::IndexType y = 0; y < m_coarserGraph->nodeCount(); ++y)
			{
				for (GraphEdge::IndexType i = coarserAdjacencyMatrix.rowBegin(y); i < coarserAdjacencyMatrix.rowEnd(y); ++i)
				{
					if (coarserAdjacencyMatrix.getX(i) == y)
					{
						coarserAdjacencyMatrix.getVs()[i] = 0;
					}
				}
			}
		}

		/*
		ここからcoarserGraphが収束した時呼ぶ関数
		*/

		// coarserGraphからノードの位置をコピーする
		template <class URBG>
		void copyFromCoarserGraph(URBG&& urbg)
		{
			Array<detail::SparseEntry<float>> perturbEntries(Arg::reserve = (nodeCount() * 2));

			for (GraphEdge::IndexType i = 0; i < nodeCount(); ++i)
			{
				const auto d = static_cast<Vector2D<float>>(RandomVec2(m_config.K * 1.e-6, std::forward<URBG>(urbg)));
				perturbEntries.emplace_back(0, i, d.x);
				perturbEntries.emplace_back(1, i, d.y);
			}

			const detail::SparseMat<float> perturb{ perturbEntries };
			m_positions = m_prolongationMatrix * m_coarserGraph->m_positions + perturb;
		}

		void prepareSimulation()
		{
			m_timeStep = m_config.initialTimeStep;
			m_progressCount = 0;
			m_stuckCount = 0;
			m_minEnergy = Largest<double>;
			m_averageEnergy = 1.0;
			m_failed = false;
			m_initialMovementSum = -1.0;
			m_currentProgressRate = getEntireProgress(0.0);
		}

		/*
		ここからupdateで毎回呼ばれる関数
		*/

		void simulateOneStep()
		{
			if (m_converged && (1 <= m_graphDepth || 0 == m_graphDepth && m_config.autoSuspend))
			{
				return;
			}

			if (nodeCount() <= 2)
			{
				m_converged = true;
				return;
			}

			const double energy0 = m_energy;
			m_energy = 0;

			reconstructQuadTree();

			const GraphEdge::IndexType threadCount = Min(static_cast<GraphEdge::IndexType>(useThreadCount()), nodeCount());
			const GraphEdge::IndexType nodeCountPerThread = nodeCount() / threadCount;

			// 並列計算での一貫性を保つために読み込み元と書き込み先は分ける
			{
				auto& oldVs = m_oldPositions.getVs();
				for (GraphEdge::IndexType nodeIndex = 0; nodeIndex < nodeCount(); ++nodeIndex)
				{
					const GraphEdge::IndexType internalIndex = m_positions.rowBegin(nodeIndex);
					const Vector2D<float> p0{ m_positions.getV(internalIndex), m_positions.getV(internalIndex + 1) };
					oldVs[internalIndex + 0] = p0.x;
					oldVs[internalIndex + 1] = p0.y;
				}
			}

			Array<std::future<TaskData>> futures;
			for (GraphEdge::IndexType i = 0; i < threadCount; ++i)
			{
				if (i + 1 == threadCount)
				{
					futures.push_back(std::async(std::launch::async, &LayoutForceDirected::updateSpringElectrical, this, i * nodeCountPerThread, nodeCount()));
				}
				else
				{
					futures.push_back(std::async(std::launch::async, &LayoutForceDirected::updateSpringElectrical, this, i * nodeCountPerThread, (i + 1) * nodeCountPerThread));
				}
			}

			double movementSum = 0.0;

			for (auto& f : futures)
			{
				const auto data = f.get();
				movementSum += data.movementSum;
				m_energy += data.energy;
			}

			if (not m_converged)
			{
				detectStuck();

				updateStepLength(energy0);
			}

			// 収束判定とprogressの更新
			{
				const double goal = m_config.K * m_config.tol;

				if (m_initialMovementSum == -1.0)
				{
					m_initialMovementSum = movementSum;
				}

				const double initialTargetScale = goal / m_initialMovementSum;
				const double currentTargetScale = goal / movementSum;

				const double initialDistFromGoal = log(initialTargetScale) / log(m_config.stepScale);
				const double currentDistFromGoal = log(currentTargetScale) / log(m_config.stepScale);
				const double localProgress = 1.0 - currentDistFromGoal / initialDistFromGoal;

				m_currentProgressRate = getEntireProgress(localProgress);
				m_converged = movementSum < m_config.K* m_config.tol;
			}
		}

		void reconstructQuadTree()
		{
			const double eps = 1.e-4;
			const RectF rect = graphBoundingRect().stretched(eps);

			m_quadTree.init(rect, m_config.maxQuadTreeDepth);

			for (auto nodeIndex : step(nodeCount()))
			{
				const GraphEdge::IndexType internalIndex = m_positions.rowBegin(nodeIndex);
				m_quadTree.add(Vec2{ m_positions.getV(internalIndex), m_positions.getV(internalIndex + 1) });
			}

			m_quadTree.update();
		}

		struct TaskData
		{
			double movementSum;
			double energy;
		};

		TaskData updateSpringElectrical(GraphEdge::IndexType startIndex, GraphEdge::IndexType endIndex)
		{
			const double p = m_config.repulsiveExponent;

			const double coeff = -m_config.C * pow(m_config.K, 1.0 + p);

			double currentEnergy = 0.0;
			double movementSum = 0.0;

			const Vec2 graphCenter = graphBoundingRect().center();

			for (GraphEdge::IndexType p0Index = startIndex; p0Index < endIndex; ++p0Index)
			{
				const GraphEdge::IndexType p0InternalIndex = m_oldPositions.rowBegin(p0Index);
				const Vec2 p0{ m_oldPositions.getV(p0InternalIndex), m_oldPositions.getV(p0InternalIndex + 1) };

				Vec2 repulsiveForce = Vec2::Zero();
				const auto fr = [=](const Vec2& xi, const Vec2& xj)
				{
					const Vec2 v = (xj - xi);
					const double lengthSq = v.lengthSq();

					return coeff * v / Math::Pow(lengthSq, p * 0.5);
				};

				if (50 < nodeCount())
				{
					// 反発力の計算: O(nlogn)
					const Vec2 f = m_quadTree.repulsiveForce(p0, p, coeff);

					if (std::isnormal(f.x) && std::isnormal(f.y))
					{
						repulsiveForce += f;
					}
				}
				else
				{
					// 反発力の計算: O(|V|^2)
					for (GraphEdge::IndexType p1Index = 0; p1Index < nodeCount(); ++p1Index)
					{
						if (p0Index == p1Index)
						{
							continue;
						}

						const GraphEdge::IndexType p1InternalIndex = m_oldPositions.rowBegin(p1Index);
						const Vec2 p1{ m_oldPositions.getV(p1InternalIndex), m_oldPositions.getV(p1InternalIndex + 1) };

						const Vec2 f = fr(p0, p1);
						if (std::isnormal(f.x) && std::isnormal(f.y))
						{
							repulsiveForce += f;
						}
					}
				}

				Vec2 attractiveForce{ 0, 0 };

				// 引力の計算
				for (GraphEdge::IndexType i = m_adjacencyMatrix.rowBegin(p0Index); i < m_adjacencyMatrix.rowEnd(p0Index); ++i)
				{
					const GraphEdge::IndexType p1Index = m_adjacencyMatrix.getX(i);

					if (p0Index == p1Index)
					{
						continue;
					}

					const GraphEdge::IndexType p1InternalIndex = m_oldPositions.rowBegin(p1Index);
					const Vec2 p1{ m_oldPositions.getV(p1InternalIndex), m_oldPositions.getV(p1InternalIndex + 1) };

					const Vec2 v = (p1 - p0);
					const Vec2 f = v.length() * v / m_config.K;
					if (std::isnormal(f.x) && std::isnormal(f.y))
					{
						attractiveForce += f;
					}
				}

				const Vec2 force = attractiveForce + repulsiveForce;

				if (1.e-7 < force.length())
				{
					Vec2 newp0 = p0 + m_timeStep * force.normalized();

					if (m_config.updateFunction)
					{
						const Vec2 p0DrawPos = toDrawPos(graphCenter, newp0);
						//const Vec2 p0DrawPos = newp0;
						const Vec2 newp0DrawPos = m_config.updateFunction(p0Index, p0, p0DrawPos);

						if (p0DrawPos != newp0DrawPos)
						{
							m_timeStep = m_config.initialTimeStep;
							newp0 = toGraphPos(graphCenter, newp0DrawPos);
							//newp0 = newp0DrawPos;
						}
					}

					movementSum += (newp0 - p0).length();

					const GraphEdge::IndexType newp0InternalIndex = m_positions.rowBegin(p0Index);
					m_positions.getVs()[newp0InternalIndex + 0] = static_cast<float>(newp0.x);
					m_positions.getVs()[newp0InternalIndex + 1] = static_cast<float>(newp0.y);

					currentEnergy += force.lengthSq();
				}
			}

			return TaskData{
				.movementSum = movementSum,
				.energy = currentEnergy,
			};
		}

		int32 coarsestDepth() const noexcept
		{
			if (not m_coarserGraph)
			{
				return m_graphDepth;
			}

			return m_coarserGraph->coarsestDepth();
		}

		// { 0.0, 1.0, gamma, gamma^2, ... gamma^maxDepth } を返す
		Array<double> getWeights() const noexcept
		{
			const double gamma = Math::Sqrt(7.0 / 4.0);
			return Array<double>::IndexedGenerate(coarsestDepth() + 2, [&](size_t i) { return i == 0 ? 0.0 : Math::Pow(gamma, i - 1); });
		}

		// depth ごとに重みを付けて progress を変換する
		double getEntireProgress(double localProgress) const
		{
			const int32 maxDepth = coarsestDepth();
			const auto weights = getWeights();

			const double weight0 = std::accumulate(weights.begin(), weights.begin() + maxDepth - m_graphDepth + 1, 0.0);
			const double weight1 = weight0 + weights[maxDepth - m_graphDepth + 1];
			const double weightSum = weights.sum();

			const double beginT = weight0 / weightSum;
			const double endT = weight1 / weightSum;

			return Math::Lerp(beginT, endT, localProgress);
		}

		void detectStuck()
		{
			if (m_energy < m_minEnergy)
			{
				m_minEnergy = m_energy;
				m_stuckCount = 0;
			}

			const int32 failStuckCount = 500;
			if (failStuckCount < m_stuckCount)
			{
				m_failed = true;
			}

			// 0.93^10 ~= 0.5
			const double history = 0.93;
			m_averageEnergy = m_averageEnergy * history + m_energy * (1.0 - history);

			const double decreaseRate = m_energy / m_averageEnergy;
			if (0.99999 < decreaseRate)
			{
				++m_stuckCount;
			}
			else
			{
				m_stuckCount = 0;
			}
		}

		// adaptive cooling scheme
		void updateStepLength(double energy0)
		{
			if (m_energy < energy0)
			{
				m_progressCount += 1;

				if (5 <= m_progressCount)
				{
					m_progressCount = 0;
					m_timeStep = m_timeStep / m_config.stepScale;
				}
			}
			else
			{
				m_progressCount = 0;
				m_timeStep = m_timeStep * m_config.stepScale;
			}
		}

		detail::SparseMat<float> m_positions;

		detail::SparseMat<float> m_oldPositions;

		detail::SparseMat<float> m_adjacencyMatrix;

		detail::SparseMat<float> m_prolongationMatrix;

		detail::SparseMat<float> m_scaleVector;

		std::unique_ptr<LayoutForceDirected> m_coarserGraph;

		// シミュレーションの実行に使う
		ForceDirectedConfig m_config;

		detail::QuadTreeVertices m_quadTree;

		double m_energy = 0.0;

		double m_timeStep = 0.0;

		int32 m_progressCount = 0;

		bool m_converged = false;


		// 失敗判定に使う
		double m_minEnergy = Largest<double>;

		double m_averageEnergy = 0.0;

		int32 m_stuckCount = 0;

		bool m_failed = false;


		// progressの計算に使う
		double m_initialMovementSum = -1.0;

		double m_currentMovementSum = 0.0;

		double m_currentProgressRate = 0.0;


		// reset時に使う
		size_t m_initRandomCount = 0;


		Array<GraphEdge::IndexType> m_originalNodeIndices;

		double m_elapsedSec = 0.0;

		int32 m_graphDepth = 0;
	};

	class LayoutCircular : public detail::GraphTransform
	{
	public:

		LayoutCircular() = default;

		explicit LayoutCircular(const ConnectedGraph& connectedGraph)
		{
			init(connectedGraph);
		}

		void init(const ConnectedGraph& connectedGraph)
		{
			m_positions.resize(connectedGraph.nodeCount());

			for (auto i : step(static_cast<GraphEdge::IndexType>(connectedGraph.nodeCount())))
			{
				m_positions[i] = Circular(Arg::r = 1.0, Arg::theta = 2_pi * i / connectedGraph.nodeCount()).fastToVec2();
			}

			m_adjacencyMatrix = detail::SparseMat<float>{ connectedGraph.edges() };
		}

		void draw(const IGraphVisualizer& visualizer) const
		{
			const Vec2 graphCenter = graphBoundingRect().center();

			for (GraphEdge::IndexType p0Index = 0; p0Index < nodeCount(); ++p0Index)
			{
				for (GraphEdge::IndexType i = m_adjacencyMatrix.rowBegin(p0Index); i < m_adjacencyMatrix.rowEnd(p0Index); ++i)
				{
					const GraphEdge::IndexType p1Index = m_adjacencyMatrix.getX(i);
					if (p0Index <= p1Index)
					{
						continue;
					}

					const Vec2& p0 = m_positions[p0Index];
					const Vec2& p1 = m_positions[p1Index];

					visualizer.drawEdge(Line{ toDrawPos(graphCenter, p0), toDrawPos(graphCenter, p1) }, p0Index, p1Index);
				}
			}

			for (GraphEdge::IndexType nodeIndex = 0; nodeIndex < nodeCount(); ++nodeIndex)
			{
				const Vec2& pos = m_positions[nodeIndex];

				visualizer.drawNode(toDrawPos(graphCenter, pos), nodeIndex);
			}
		}

		[[nodiscard]]
		GraphEdge::IndexType nodeCount() const noexcept
		{
			return static_cast<GraphEdge::IndexType>(m_positions.size());
		}

		[[nodiscard]]
		GraphEdge::IndexType edgeCount() const noexcept
		{
			return static_cast<GraphEdge::IndexType>(m_adjacencyMatrix.getXs().size());
		}

		[[nodiscard]]
		RectF boundingRect() const
		{
			const auto graphRect = graphBoundingRect();
			const Vec2 tl = toDrawPos(graphRect.center(), graphRect.tl());
			const Vec2 br = toDrawPos(graphRect.center(), graphRect.br());
			return RectF{ tl, br - tl };
		}

		[[nodiscard]]
		Vec2 centroid() const
		{
			const auto graphRect = graphBoundingRect();
			const auto pos = graphCentroid();
			return toDrawPos(graphRect.center(), pos);
		}

	private:

		[[nodiscard]]
		RectF graphBoundingRect() const override
		{
			return Geometry2D::BoundingRect(m_positions);
		}

		[[nodiscard]]
		Vec2 graphCentroid() const
		{
			if (nodeCount() == 0)
			{
				return Vec2{};
			}

			return m_positions.sum() / static_cast<double>(m_positions.size());
		}

		Array<Vec2> m_positions;

		detail::SparseMat<float> m_adjacencyMatrix;
	};
}
