# OpenSiv3D_GraphDrawSample

OpenSiv3D でグラフ描画を行うサンプルとチュートリアルです。

# サンプル

<table align="center">
    <tr>
        <td align="center"><a href="https://github.com/agehama/OpenSiv3D_GraphDrawSample/blob/main/example/example1_TextGraph.cpp"><b>TextGraph</b></a></td>
	<td align="center"><a href="https://github.com/agehama/OpenSiv3D_GraphDrawSample/blob/main/example/example2_MultipleGraphs.cpp"><b>MultipleGraphs</b></a></td>
    </tr>
    <tr>
        <td align="center"><img alt="example1" src="https://user-images.githubusercontent.com/4939010/125158652-5ebd7b00-e1ad-11eb-99d0-a1b0b0c5cd3e.png" width="600px"></td>
        <td align="center"><img alt="example2" src="https://user-images.githubusercontent.com/4939010/125158654-611fd500-e1ad-11eb-9e58-05c9d2cd3246.png" width="600px"></td>
    </tr>
    <tr>
        <td align="center"><a href="https://github.com/agehama/OpenSiv3D_GraphDrawSample/blob/main/example/example3_JSONViewer.cpp"><b>JSONViewer</b></a></td>
	<td align="center"><a href="https://github.com/agehama/OpenSiv3D_GraphDrawSample/blob/main/example/example4_PathSearch.cpp"><b>PathSearch</b></a></td>
    </tr>
    <tr>
        <td align="center"><img alt="example3" src="https://user-images.githubusercontent.com/4939010/125158658-62510200-e1ad-11eb-9657-61ac78c0a1c1.png" width="600px"></td>
        <td align="center"><img alt="example4" src="https://user-images.githubusercontent.com/4939010/125159539-025d5a00-e1b3-11eb-804c-d717ebf368a1.png" width="600px"></td>
    </tr>
</table>

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
	const GraphLoader loader(U"example/primitives.mtx");

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
	const GraphLoader loader{ U"example/simpleGraph.txt" };

	auto layout = LayoutForceDirected{ loader[0], ForceDirectedConfig{ .startImmediately = StartImmediately::Yes } };

	RectF rect = Scene::Rect().stretched(-100);

	const BasicGraphVisualizer visualizer;

	while (System::Update())
	{
		rect.drawFrame(2.0);

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
	RectF rect = Scene::Rect().stretched(-100);

-	const BasicGraphVisualizer visualizer;
+	Scene::SetBackground(Color(U"#f7f1cf"));
+	BasicGraphVisualizer visualizer{ 15, 5, Color(U"#7adb6b"), Color(U"#e5da9a") };

	while (System::Update())
```

また、`layout.setDrawArea()` は渡された `Rect` にノードの中心座標を揃えるため、(1) のプログラムは描画したときに外側のノードがはみ出ています。
これを描画範囲にぴったり収めるにはノードの半径分縮めた `Rect` を `layout.setDrawArea()` に渡すようにします。
```diff
		}

-		layout.setDrawArea(rect);
+		layout.setDrawArea(rect.stretched(-visualizer.m_nodeRadius));

		layout.draw(visualizer);
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
	Scene::SetBackground(Color(U"#f7f1cf"));
-	BasicGraphVisualizer visualizer{ 15, 5, Color(U"#7adb6b"), Color(U"#e5da9a") };
+	LabelGraphVisualizer visualizer{ Font{16, Typeface::Heavy }, Color(U"#f7f1cf"), 15, 5, Color(U"#7adb6b"), Color(U"#e5da9a") };

	while (System::Update())
```

<p align="center">
  <img alt="tutorial_2_3" src="https://user-images.githubusercontent.com/4939010/121796413-bc64b300-cc53-11eb-894f-6ccdd7f4a53d.png" width="60%">
</p>

### (4) レイアウトを固定する

`LayoutForceDirected` は乱数を使うため実行するたびに異なるレイアウトに収束します。
予め乱数のシードを固定することで、同じレイアウトを再現することが可能です。

```diff
	const GraphLoader loader{ U"example/simpleGraph.txt" };

+	GetDefaultRNG().seed(0); // シード値を0に設定

	auto layout = LayoutForceDirected{ loader[0], ForceDirectedConfig{.startImmediately = StartImmediately::Yes } };
```

### (5) 回転する

`Transformer2D` を作って描画範囲ごと回転したりスケールをかけたりすることができます。

```diff
+	double angle = 30_deg;

	while (System::Update())
	{
+		// マウスホイールで回転する
+		angle += Mouse::Wheel() * 0.1;
+
+		const auto mat = Mat3x2::Rotate(angle, Scene::Center());
+		const Transformer2D t(mat, TransformCursor::Yes);

		rect.drawFrame(2.0);
```

<p align="center">
  <img alt="tutorial_2_5_1" src="https://user-images.githubusercontent.com/4939010/122674134-c156d300-d20e-11eb-97bb-d75eb6cfc8f0.gif" width="60%">
</p>

ここで `layout.setDrawArea()` の第二引数に `Mat3x2` を渡せば、描画範囲を固定したままトランスフォームをかけることができます。

```diff
	while (System::Update())
	{
		// マウスホイールで回転する
		angle += Mouse::Wheel() * 0.1;

+		rect.drawFrame(2.0);
+
		const auto mat = Mat3x2::Rotate(angle, Scene::Center());
-		const Transformer2D t(mat, TransformCursor::Yes);
+		const Transformer2D t(mat);

-		rect.drawFrame(2.0);
```

```diff
		}

-		layout.setDrawArea(rect.stretched(-visualizer.m_nodeRadius));
+		layout.setDrawArea(rect.stretched(-visualizer.m_nodeRadius), mat);

		layout.draw(visualizer);
	}
```

<p align="center">
  <img alt="tutorial_2_5_2" src="https://user-images.githubusercontent.com/4939010/122674888-f1ec3c00-d211-11eb-86e0-6f269e111552.gif" width="60%">
</p>

### チュートリアル2 ソースコード全体

```cpp
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

	virtual void drawEdge(const Line& line, GraphEdge::IndexType, GraphEdge::IndexType) const override
	{
		line.draw(LineStyle::RoundDot, m_edgeThickness, m_edgeColor);
	}

	Font m_labelFont;

	ColorF m_labelColor;
};

void Main()
{
	const GraphLoader loader{ U"example/simpleGraph.txt" };

	GetDefaultRNG().seed(0); // シード値を0に設定

	auto layout = LayoutForceDirected{ loader[0], ForceDirectedConfig{.startImmediately = StartImmediately::Yes } };

	RectF rect = Scene::Rect().stretched(-100);

	Scene::SetBackground(Color(U"#f7f1cf"));
	LabelGraphVisualizer visualizer{ Font{16, Typeface::Heavy }, Color(U"#f7f1cf"), 15, 5, Color(U"#7adb6b"), Color(U"#e5da9a") };

	double angle = 30_deg;

	while (System::Update())
	{
		// マウスホイールで回転する
		angle += Mouse::Wheel() * 0.1;

		rect.drawFrame(2.0);

		const auto mat = Mat3x2::Rotate(angle, Scene::Center());
		const Transformer2D t(mat);

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

		layout.setDrawArea(rect.stretched(-visualizer.m_nodeRadius), mat);

		layout.draw(visualizer);
	}
}
```

## 3 応用編 インタラクティブな描画

ForceDirected レイアウトを使ってグラフの配置をインタラクティブに編集するアプリケーションを作ってみます。

### (1) ループでグラフをレイアウトする

これまでは全て初期化時にレイアウトの計算を行っていました。
しかし、規模の大きいグラフ（ノード数が10000以上）になるとレイアウトの計算が完了するまでに数十秒かかることもあります。

このような場合、初期化時に計算を行わずにループの中で `layout.update()` を呼ぶことで、描画を更新しながらレイアウトの計算を行うことができます。
`layout.update()` の引数には、計算に使う時間をミリ秒で指定します。
他に重い処理がないプログラムであれば、`16` ミリ秒としておけば 60FPS を維持しながら計算を進めます。

また、これまでは `layout.setDrawArea()` は最初に一度呼んだきりでしたが、レイアウトが更新されるたびに座標が変わるので呼びなおす必要があります。

```cpp
void Main()
{
	const GraphLoader loader{ U"example/sierpinski.txt" };

	const double nodeRadius = 7;
	BasicGraphVisualizer visualizer{ nodeRadius };

	GetDefaultRNG().seed(0);

	LayoutForceDirected layout{ loader[0], ForceDirectedConfig{} };

	while (System::Update())
	{
		layout.update(16);

		layout.setDrawArea(Scene::Rect().stretched(-50));

		layout.draw(visualizer);
	}
}
```

<p align="center">
  <img alt="tutorial_3_1" src="https://user-images.githubusercontent.com/4939010/122638885-bda14e80-d131-11eb-9c53-bcfb083209db.png" width="60%">
</p>

### (2) ノードのマウスクリックを実装する

ノードの現在位置は `layout.activeNodePositions()` で取得することができます。

これを使ってクリックされたノードのインデックスを表示する機能を追加します。

```diff
void Main()
{
+	const Font font(16, Typeface::Heavy);
+
+	Optional<GraphEdge::IndexType> clickedNode;

	const GraphLoader loader{ U"example/sierpinski.txt" };
```

```diff
	while (System::Update())
	{
		layout.update(16);

		layout.setDrawArea(Scene::Rect().stretched(-50));

		layout.draw(visualizer);

+		for (auto& [nodeIndex, nodePos] : layout.activeNodePositions())
+		{
+			const auto nodeCircle = nodePos.asCircle(nodeRadius);
+
+			if (nodeCircle.leftClicked())
+			{
+				clickedNode = nodeIndex;
+			}
+
+			if (clickedNode == nodeIndex)
+			{
+				nodeCircle.draw(Palette::Orange);
+
+				const auto labelPos = nodePos + Circular{ 20, 30_deg };
+				font(nodeIndex).draw(labelPos + Vec2{ 1, 1 }, Palette::Black);
+				font(nodeIndex).draw(labelPos, Palette::Orange);
+			}
+		}
	}
```

<p align="center">
  <img alt="tutorial_3_2" src="https://user-images.githubusercontent.com/4939010/122639035-a020b480-d132-11eb-9194-ca688f32726c.gif" width="60%">
</p>


### (3) マウスでドラッグしてノードを動かす

ノードのクリック判定が取れるようになったので、次はクリックしたノードの座標をカーソル位置に移動するようにします。

まず `config.autoSuspend` を `false` にしてレイアウトが完了しても座標更新を続けるようにします。
そして `config.updateFunction` にはレイアウト計算でそれぞれのノードに対して呼ばれる座標更新関数を設定することができます。
これにクリック中のノードの座標をカーソル位置に移動する処理を加えましょう。

```diff
	GetDefaultRNG().seed(0);

-	LayoutForceDirected layout{ loader[0], ForceDirectedConfig{} };
+	ForceDirectedConfig config
+	{
+		.autoSuspend = false,
+		.initialTimeStep = 0.01, // クリック時の見た目のぶれを小さくするため
+	};
+
+	config.updateFunction = [&](GraphEdge::IndexType nodeIndex, const Vec2& /*oldPos*/, const Vec2& newPos)
+	{
+		if (clickedNode && clickedNode.value() == nodeIndex)
+		{
+			return Cursor::PosF();
+		}
+
+		return newPos;
+	};
+
+	LayoutForceDirected layout{ loader[0], config };
```

あとはマウスを離したときに `clickedNode` をリセットする処理を入れればドラッグ移動ができるようになります。

ただし、これだけだとドラッグしながら描画範囲から出た時に `layout.setDrawArea()` で全体を縮小する処理が連続して走ってしまうため、ドラッグ中は `layout.setDrawArea()` が呼ばれないように変更します。

```diff
	while (System::Update())
	{
		layout.update(16);
		
-		layout.setDrawArea(Scene::Rect().stretched(-50));
+		if (!MouseL.pressed())
+		{
+			clickedNode = none;
+
+			layout.setDrawArea(Scene::Rect().stretched(-50));
+		}

		layout.draw(visualizer);
```

<p align="center">
  <img alt="tutorial_3_3" src="https://user-images.githubusercontent.com/4939010/122639262-ee828300-d133-11eb-8723-736e7d22b1a3.gif" width="60%">
</p>


### チュートリアル3 ソースコード全体

```cpp
#include <Siv3D.hpp> // OpenSiv3D v0.6

#include "include/GraphDrawing.hpp"

void Main()
{
	const Font font(16, Typeface::Heavy);

	Optional<GraphEdge::IndexType> clickedNode;

	const GraphLoader loader{ U"example/sierpinski.txt" };

	const double nodeRadius = 7;
	BasicGraphVisualizer visualizer{ nodeRadius };

	GetDefaultRNG().seed(0);

	ForceDirectedConfig config
	{
		.autoSuspend = false,
		.initialTimeStep = 0.01, // クリック時の見た目のぶれを小さくするため
	};

	config.updateFunction = [&](GraphEdge::IndexType nodeIndex, const Vec2& /*oldPos*/, const Vec2& newPos)
	{
		if (clickedNode && clickedNode.value() == nodeIndex)
		{
			return Cursor::PosF();
		}

		return newPos;
	};

	LayoutForceDirected layout{ loader[0], config };

	while (System::Update())
	{
		layout.update(16);

		if (!MouseL.pressed())
		{
			clickedNode = none;

			layout.setDrawArea(Scene::Rect().stretched(-50));
		}

		layout.draw(visualizer);

		for (auto& [nodeIndex, nodePos] : layout.activeNodePositions())
		{
			const auto nodeCircle = nodePos.asCircle(nodeRadius);

			if (nodeCircle.leftClicked())
			{
				clickedNode = nodeIndex;
			}

			if (clickedNode == nodeIndex)
			{
				nodeCircle.draw(Palette::Orange);

				const auto labelPos = nodePos + Circular{ 20, 30_deg };
				font(nodeIndex).draw(labelPos + Vec2{ 1, 1 }, Palette::Black);
				font(nodeIndex).draw(labelPos, Palette::Orange);
			}
		}
	}
}
```
