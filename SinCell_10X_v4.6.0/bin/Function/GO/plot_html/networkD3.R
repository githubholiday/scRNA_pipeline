#! /annoroad/share/software/install/R-3.2.2/bin/Rscript
library(networkD3)
library(magrittr)
library(htmlwidgets)
args=commandArgs(T)
node <- args[1]
link <- args[2]
outfile <- args[3]

links<-read.table(link,sep='\t',header=TRUE)
nodes<-read.table(node,sep='\t',header=TRUE)
fun <- forceNetwork(Links = links,#线性质数据框
		Nodes = nodes,#节点性质数据框
		#数据库的标题 
		Source = "source",#连线的源变量
		Target = "target",#连线的目标变量
		Value = "value",#连线的粗细值
		NodeID = "name",#节点名称
		Group = "group",#节点的分组
		Nodesize = "size" ,#节点大小，节点数据框中
		#设置图框宽度和高度
	#	width = 850,
	#	height = 480,
		###美化部分
          #   fontFamily="宋体",#字体设置如"华文行楷" 等
		fontSize = 24, #节点文本标签的数字字体大小（以像素为单位）。
            # linkColour="black",#连线颜色,black,red,blue, 
		colourScale=JS('d3.scaleOrdinal().domain(["up", "down"]).range(["green","blue"])'),
		#colourScale ,linkWidth,#节点颜色,red，蓝色blue,cyan,yellow等
		charge = -100,#数值表示节点排斥强度（负值）或吸引力（正值）  
		opacity = 1.0, #所有节点初始透明度
		legend=T,#显示节点分组的颜色标签
		arrows=F,#是否带方向
	#	bounded=T,#是否启用限制图像的边框
		opacityNoHover=0,
		zoom = T)
htmlwidgets::onRender(fun,'function(el,x) { 
    var link = d3.selectAll(".link")
    var node = d3.selectAll(".node")

    var options = { opacity: 1,
                    clickTextSize: 10,
                    opacityNoHover: 0.1,
                    radiusCalculation: "Math.sqrt(d.nodesize)+6"
                  }

    var unfocusDivisor = 4;

    var links = HTMLWidgets.dataframeToD3(x.links);
    var linkedByIndex = {};

    links.forEach(function(d) {
      linkedByIndex[d.source + "," + d.target] = 1;
      linkedByIndex[d.target + "," + d.source] = 1;
    });

    function neighboring(a, b) {
      return linkedByIndex[a.index + "," + b.index];
    }
function nodeSize(d) {
            if(options.nodesize){
                    return eval(options.radiusCalculation);
            }else{
                    return 6}

    }

    function mouseover(d) {
      var unfocusDivisor = 4;

      link.transition().duration(200)
        .style("opacity", function(l) { return d != l.source && d != l.target ? +options.opacity / unfocusDivisor : +options.opacity });

      node.transition().duration(200)        .style("opacity", function(o) { return d.index == o.index || neighboring(d, o) ? +options.opacity : +options.opacity / unfocusDivisor; });

      d3.select(this).select("circle").transition()
        .duration(750)
        .attr("r", function(d){return nodeSize(d)+5;});

      node.select("text").transition()
.duration(750)
        .attr("x", 13)
        .style("stroke-width", ".5px")
        .style("font", 24 + "px ")
        .style("opacity", function(o) { return d.index == o.index || neighboring(d, o) ? 1 : 0; });
    }

    function mouseout() {
      node.style("opacity", +options.opacity);
      link.style("opacity", +options.opacity);

      d3.select(this).select("circle").transition()
        .duration(750)
        .attr("r", function(d){return nodeSize(d);});
      node.select("text").transition()
        .duration(1250)
        .attr("x", 0)
        .style("font", options.fontSize + "px ")
        .style("opacity", 0);
    }

    d3.selectAll(".node").on("mouseover", mouseover).on("mouseout", mouseout);
}')%>%
saveNetwork(file = outfile,selfcontained=F)

