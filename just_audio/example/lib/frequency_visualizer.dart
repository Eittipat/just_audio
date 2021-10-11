import 'dart:core';

import 'package:flutter/material.dart';
import 'package:flutter/services.dart';
import 'package:flutter/widgets.dart';
import 'package:just_audio/just_audio.dart';

class FrequencyVisualizerWidget extends StatefulWidget {
  final VisualizerFftCapture capture;

  FrequencyVisualizerWidget(this.capture);
  @override
  FrequencyVisualizerWidgetState createState() => FrequencyVisualizerWidgetState();
}

class FrequencyVisualizerWidgetState extends State<FrequencyVisualizerWidget> {
  @override
  Widget build(BuildContext context) {
    return ClipRect(
      child: CustomPaint(
        painter: FrequencyVisualizerPainter(
          capture: widget.capture,
          color: Colors.blueAccent,
        ),
      ),
    );
  }
}

// https://github.com/iamSahdeep/FlutterVisualizers/blob/master/lib/Visualizers/BarVisualizer.dart
class FrequencyVisualizerPainter extends CustomPainter {
  final VisualizerFftCapture capture;
  final Color color;
  final Paint wavePaint;
  final int density;
  final int gap;

  FrequencyVisualizerPainter({
    required this.capture,
    required this.color,
    this.density = 100,
    this.gap = 2,
  }) : wavePaint = Paint()
          ..color = color.withOpacity(1.0)
          ..style = PaintingStyle.fill;

  @override
  void paint(Canvas canvas, Size size) {
    final width = size.width;
    final height = size.height;
    final buffer = capture.data.buffer;
    final bytes = ByteData.view(buffer);

    double barWidth = width / density;
    double div = capture.data.length / density;
    wavePaint.strokeWidth = barWidth - gap;
    for (int i = 0; i < density; i++) {
      int bytePosition = (i * div).ceil();
      int value = bytes.getInt8(bytePosition);
      double top = (height / 2 - (((value)).abs()));
      double barX = (i * barWidth) + (barWidth / 2);
      if (top > height) top = top - height;
      canvas.drawLine(Offset(barX, height / 2), Offset(barX, top), wavePaint);
    }
  }

  @override
  bool shouldRepaint(CustomPainter oldDelegate) {
    return true;
  }
}
