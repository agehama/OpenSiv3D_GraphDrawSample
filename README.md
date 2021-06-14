# OpenSiv3D_GraphDrawSample

OpenSiv3D でグラフ描画を行うサンプルとチュートリアルです。

グラフ描画とは、グラフのデータ構造（何番目のノード同士がつながっているかという接続関係）をもとに、適切なノードの座標を計算して可視化する手法のことです。

# チュートリアル

## 1 読み込みと初期化

### (1) Circular レイアウトで描画する

まず、描画するグラフのエッジの配列を引数にして `ConnectedGraph` クラスを作ります（`ConnectedGraph` はすべてのノードがエッジでつながれた`連結グラフ`を表すデータ構造です。）

次に `ConnectedGraph` を `LayoutCircular` クラスに渡して初期化します。
ここで円形にノードを配置する座標が計算されます。

その後 `layout.draw()` を呼んで画面にグラフを描画します。
このとき、`LayoutCircular` はグラフの配置情報しか持たないので、引数に `BasicGraphVisualizer` クラスを与えて色や大きさなど描画の方法を指定します。

```cpp
void Main()
{
	const ConnectedGraph graph = { {
		{0, 1},
		{2, 1},
		{1, 3},
		{3, 4},
		{3, 5},
		{4, 6},
		{5, 6},
	} };

	const LayoutCircular layout{ graph };

	while (System::Update())
	{
		layout.draw(BasicGraphVisualizer{});
	}
}
```

<p align="center">
  <img alt="tutorial_1_1" src="https://user-images.githubusercontent.com/4939010/121796481-4b71cb00-cc54-11eb-8295-1e2fbbfd1285.png" width="60%">
</p>

### (2) ForceDirected レイアウトで描画する

今度は `LayoutForceDirected` クラスに `ConnectedGraph` を渡して  ForceDirected レイアウトを行います。

`LayoutForceDirected` のレイアウト計算は複雑なグラフに対して時間がかかるため、Circular レイアウトと異なり通常は `layout.update()` を呼んだタイミングでのみ行われます。

ここでは、例として簡単なグラフを扱うので、設定に `.startImmediately = StartImmediately::Yes` を指定してレイアウト計算をその場で実行します（複雑なグラフをループで少しずつ計算する方法は [3 応用編 インタラクティブな描画](#3-応用編-インタラクティブな描画) を参照してください。）

```cpp
void Main()
{
	const ConnectedGraph graph = { {
		{0, 1},
		{2, 1},
		{1, 3},
		{3, 4},
		{3, 5},
		{4, 6},
		{5, 6},
	} };

	const LayoutForceDirected layout{ graph, ForceDirectedConfig{ .startImmediately = StartImmediately::Yes } };

	while (System::Update())
	{
		layout.draw(BasicGraphVisualizer{});
	}
}
```

<p align="center">
  <img alt="tutorial_1_2" src="https://user-images.githubusercontent.com/4939010/121796483-50367f00-cc54-11eb-912d-d14b16025eb7.png" width="60%">
</p>

### (3) ファイルからグラフを読み込む

ファイルからグラフを読み込むには `GraphLoader` クラスを使います。

`GraphLoader` は入力に以下の形式をサポートします。
- Matrix Market Exchange Formats 形式 (.mtx)
- エッジリスト (.txt)
- `Array<GraphEdge>`

`GraphLoader` で読み込んだグラフは連結成分ごとに分解されて、それぞれの `ConnectedGraph` には添え字アクセスすることができます。

```cpp
void Main()
{
	const GraphLoader loader(U"primitives.mtx");

	const ForceDirectedConfig config{
		.startImmediately = StartImmediately::Yes,
	};

	int32 index = 0;
	LayoutForceDirected layout{ loader[index], config };

	const Font font{ 24 };

	while (System::Update())
	{
		if (KeySpace.down())
		{
			index = (index + 1) % loader.size();

			layout = LayoutForceDirected{ loader[index], config };
		}

		layout.draw(BasicGraphVisualizer{});

		font(U"グラフ", index + 1, U"/", loader.size(), U"（Space キーでグラフを切り替える）").draw(0, 0, Palette::Yellow);
	}
}
```

<p align="center">
  <img alt="tutorial_1_3" src="https://user-images.githubusercontent.com/4939010/121870409-fca06000-cd3d-11eb-912a-11eca0e95348.gif" width="60%">
</p>

## 2 配置と描画

### (1) Rect で指定した範囲に描画する

`layout.setDrawArea()` で描画する範囲を指定することができます。

例として `Rect` の端を掴んで描画範囲を動かせるプログラムを作ってみます。

```cpp
void Main()
{
	const GraphLoader loader{ U"example1.txt" };

	auto layout = LayoutForceDirected{ loader[0], ForceDirectedConfig{ .startImmediately = StartImmediately::Yes } };

	RectF rect = Scene::Rect().stretched(-100);

	const BasicGraphVisualizer visualizer;

	while (System::Update())
	{
		const Circle cursorCircle{ Cursor::Pos(), 30.0 };

		const bool mouseOverLeft = rect.left().intersects(cursorCircle);
		const bool mouseOverRight = rect.right().intersects(cursorCircle);
		const bool mouseOverTop = rect.top().intersects(cursorCircle);
		const bool mouseOverBottom = rect.bottom().intersects(cursorCircle);

		Cursor::SetDefaultStyle(CursorStyle::Default);

		if (mouseOverLeft || mouseOverRight)
		{
			Cursor::SetDefaultStyle(CursorStyle::ResizeLeftRight);
		}
		else if (mouseOverTop || mouseOverBottom)
		{
			Cursor::SetDefaultStyle(CursorStyle::ResizeUpDown);
		}

		if (MouseL.pressed())
		{
			if (mouseOverLeft)
			{
				rect = RectF(Arg::bottomRight = rect.br(), rect.br().x - Cursor::Pos().x, rect.h);
			}
			else if (mouseOverRight)
			{
				rect = RectF(Arg::topLeft = rect.tl(), Cursor::Pos().x - rect.tl().x, rect.h);
			}
			else if (mouseOverTop)
			{
				rect = RectF(Arg::bottomRight = rect.br(), rect.w, rect.br().y - Cursor::Pos().y);
			}
			else if (mouseOverBottom)
			{
				rect = RectF(Arg::topLeft = rect.tl(), rect.w, Cursor::Pos().y - rect.tl().y);
			}
		}

		layout.setDrawArea(rect);

		rect.drawFrame(2.0);

		layout.draw(visualizer);
	}
}
```

<p align="center">
  <img alt="tutorial_2_1" src="https://user-images.githubusercontent.com/4939010/121796403-a7881f80-cc53-11eb-8046-34e89ab52821.gif" width="60%">
</p>



### (2) 色を変える

`BasicGraphVisualizer` の引数にノードの半径、エッジの太さ、ノードの色、エッジの色を指定することができます。

```diff
-	const BasicGraphVisualizer visualizer;
+	Scene::SetBackground(Color(U"#f7f1cf"));
+	BasicGraphVisualizer visualizer{ 15, 5, Color(U"#7adb6b"), Color(U"#e5da9a") };
```

また、`layout.setDrawArea()` は渡された `Rect` にノードの中心座標を揃えるため、(1) のプログラムは描画したときに外側のノードがはみ出ています。
これを描画範囲にぴったり収めるにはノードの半径分縮めた `Rect` を `layout.setDrawArea()` に渡すようにします。
```diff
-		layout.setDrawArea(rect);
+		layout.setDrawArea(rect.stretched(-visualizer.m_nodeRadius));
```

<p align="center">
  <img alt="tutorial_2_2" src="https://user-images.githubusercontent.com/4939010/121796412-b7076880-cc53-11eb-8acb-51c56eece551.png" width="60%">
</p>

### (3) ラベルを付ける

`BasicGraphVisualizer` クラスを継承して描画関数をカスタマイズすることができます。

`drawNode()` をオーバーライドしてノードの描画をラベル付きにする `LabelGraphVisualizer` クラスを作ってみましょう。

ついでに `drawEdge()` も書き換えてエッジの描画スタイルを点線に変更してみます。

```cpp
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

	virtual void drawEdge(const Line& line, GraphEdge::IndexType, GraphEdge::IndexType) const override
	{
		line.draw(LineStyle::RoundDot, m_edgeThickness, m_edgeColor);
	}

	Font m_labelFont;

	ColorF m_labelColor;
};
```

そして `visualizer` を上で定義した `LabelGraphVisualizer` に置き換えます。

```diff
-	BasicGraphVisualizer visualizer{ 15, 5, Color(U"#7adb6b"), Color(U"#e5da9a") };
+	LabelGraphVisualizer visualizer{ Font{16, Typeface::Heavy }, Color(U"#f7f1cf"), 15, 5, Color(U"#7adb6b"), Color(U"#e5da9a") };
```

<p align="center">
  <img alt="tutorial_2_3" src="https://user-images.githubusercontent.com/4939010/121796413-bc64b300-cc53-11eb-894f-6ccdd7f4a53d.png" width="60%">
</p>

### (4) 回転する

`Transformer2D` を作って描画範囲ごと回転したりスケールをかけたりすることができます。

```diff

+	double angle = 30_deg;

	while (System::Update())
	{
+		if (KeyLeft.pressed())
+		{
+			angle -= 1_deg;
+		}
+		else if (KeyRight.pressed())
+		{
+			angle += 1_deg;
+		}

+		const Transformer2D t(Mat3x2::Rotate(angle, Scene::Center()), TransformCursor::Yes);

		const Circle cursorCircle{ Cursor::Pos(), 30.0 };
```

<p align="center">
  <img alt="tutorial_2_4" src="https://user-images.githubusercontent.com/4939010/121796414-bff83a00-cc53-11eb-8c99-e228225006aa.png" width="60%">
</p>

## 3 応用編 インタラクティブな描画
