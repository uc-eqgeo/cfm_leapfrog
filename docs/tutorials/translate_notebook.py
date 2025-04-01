import json

# Load the original notebook
with open('define_connected_no_smoothing_japanese.ipynb', 'r', encoding='utf-8') as f:
    notebook = json.load(f)

# Translated comments
translated_comments = [
    [
        "## スムージングなしで接続された断層システムを定義する\n",
        "このノートブックは、シェープファイル内のどの断層セグメントが大きな断層システムに接続されるべきかを識別する方法を示しています。これは[標準的なワークフロー](https://uc-eqgeo.github.io/cfm_leapfrog/tutorials/define_connected.html)と同じですが、いくつかのコマンドに追加の引数があります。**入力GISデータに既に密接に配置された頂点がある場合を除き、スムージングを使用することをお勧めします**。この代替ワークフローは、1kmの頂点間隔を持つNZ CFMのバージョンに対応するために実装されました。\n"
    ],
    [
        "## 断層の読み込み\n",
        "まず、断層のGIS表現を読み込む必要があります。この例では、ニュージーランドコミュニティ断層モデル（Seebeck et al., 2022）の断層のサブセットを使用します。\n"
    ],
    [
        "## セグメント間の接続を見つける\n",
        "データが読み込まれたら、距離許容値を設定する必要があります。この許容値は、接続としてカウントされる2つの断層トレース間の最小水平距離です。\n"
    ],
    [
        "次のセルは、指定された距離許容値内にあるセグメントトレースを見つけるためにpythonモジュールnetworkxを使用します。\n"
    ],
    [
        "これらの接続を手動で編集およびレビューするために書き出す必要があります。ファイルはこのJupyterノートブックと同じディレクトリに書き出されます。接頭辞はあなたが指定し、接尾辞は\"_suggested.csv\"です。\n"
    ],
    [
        "これにより、次のようなCSVファイルが作成されます：\n",
        "![unedited_connected](./tutorial_images/connected_unedited.png)\n",
        "結合された断層システムの名前は最初の列にあり、接続されたシステムを構成する断層の名前は後続の列にあります。\n"
    ],
    [
        "## 手動編集の実施と組み込み\n",
        "自動生成された断層システムの提案には（設計上）分割する必要がある過剰に接続された断層システムが含まれます。この段階では、これらのネットワークを小さな断層システムに分割する最良の方法は、CSVファイルを手動で編集することです。以下の例は、Hope Faultシステムを表す新しい行をCSVに追加したものです -- Hope Faultは自動生成された接続CSVでAlpineおよびKekerengu-Needles断層システムとグループ化されています。**この新しいファイルを保存する際には、上書きしないように異なる名前で保存してください！**\n",
        "![edited_connected](./tutorial_images/connected_edited.png)\n",
        "必要な編集を行ったら、新しいCSVを読み込んで自動生成された接続断層システムを上書きします：\n"
    ],
    [
        "## 切断階層の定義\n",
        "断層システムを定義したら（後でメッシュを作成するために使用します）、他の断層に対してどの断層が終了するかを指定する必要があります。例えば、Hope Faultの西端が深さでAlpine Faultによって切断される可能性が高いようです。この複雑なメッシュ切断は、leapfrogなどの専用ソフトウェアを使用して最も効果的に行われますが、自動的に切断を行うためには、最初に*切断階層*を指定するのが最善です。\n"
    ],
    [
        "## すべり速度に基づく階層の提案\n",
        "階層の最初のパスは、断層のすべり速度に基づいて純粋に生成できます。既に読み込んだ断層データにすべり速度が関連付けられていると仮定すると、この最初のパスは簡単に行えます：\n",
        "この操作は、モデル内の断層（または断層システム）をすべり速度の降順に並べ替えるだけです。異なるすべり速度を持つセグメントを持つ接続された断層システムの場合、その断層システム内の任意のセグメントの最大すべり速度が使用され、切断階層に断層システムが配置されます。\n"
    ],
    [
        "## 切断階層の編集\n",
        "次に、ファイル内の行の順序を切り替えることでこの階層を編集できます。交差する断層/システムのペアについては、ファイルの下部に近い断層が上部に近い断層に対して終了します。\n",
        "編集が望ましい状況の例を以下に示します。Jordan-Kekerengu-Needles Fault Systemの最大すべり速度（23 mm/yr）はHope Faultの対応する最大値（15.8 mm/yr）よりも速いですが、Jordan FaultがHope Faultに対して終了する断層モデルを作成したいと考えています。この終了を実現するために、CSVファイル内で`Hope combined`を`Jordan - Kek - Needles combined`の上に移動します。同様の理由で、`Hanmer` Faultを`Hope Hanmer NW`の上に移動します。\n",
        "![Adjust hierarchy](./tutorial_images/reorder_hierarchy.png)\n",
        "次のようにしてこの新しい階層を読み込みます：\n"
    ],
    [
        "## メッシュ作成のためのシェープファイルの作成\n",
        "メッシュ作成前の最終ステップは、断層の三角メッシュ表現を作成するためにメッシュ作成ソフトウェアと組み合わせることができるファイルの作成です。\n",
        "これらの三角形の表面を複数のソフトウェアパッケージ（例えば、MOVE 3D）で構築することは可能ですが、以下の議論はLeapfrog Geoソフトウェアの使用を前提としています。\n",
        "## シェープファイルを保持するディレクトリの作成\n",
        "組織上の理由から、異なるシェープファイルを異なるディレクトリに配置することが役立ちます\n",
        "## シェープファイルの書き出し\n"
    ]
]

# Replace the original comments with the translated comments
cell_i = 0
for cell in notebook['cells']:
    if cell['cell_type'] == 'markdown':
        cell['source'] = translated_comments[cell_i]
        cell_i += 1

# Save the translated notebook
with open('define_connected_no_smoothing_japanese_translated.ipynb', 'w', encoding='utf-8') as f:
    json.dump(notebook, f, ensure_ascii=False, indent=2)