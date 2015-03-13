#ifndef GLWIDGET_H
#define GLWIDGET_H

#include "Common.h"
#include "Lattice.h"
#include "Mouse.h"
#include <QtGui>
#include <QGLWidget>
#include <QBrush>
#include <QFont>
#include <QPen>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QVector2D>

class QPainter;
class QPaintEvent;
class QWidget;

enum drawModes {
	testParticles = 1, vectorField, vectorFieldRenormalized, fourierField
};

class GLWidget : public QGLWidget {
	Q_OBJECT

	public:
		GLWidget(QWidget *parent);
		void mouseMoveEvent(QMouseEvent *event);
		void keyPressEvent(QKeyEvent *event);

	public slots:
		void animate();

	protected:
		void paintEvent(QPaintEvent *event);

	private:
		void paint(QPainter *painter, QPaintEvent *event);
		void paintParticles(QPainter *painter);
		void paintObstacles(QPainter *painter);
		void paintVectorField(QPainter *painter);
		void paintFourierField(QPainter *painter);
		void drawArrow(QPainter *painter, QPoint p0, QPointF scale, double angle);
		
		QBrush background;
		QPen arrowPen;
		QPen particlePen;
		QColor arrowColor;
		QColor particleColor;
		QColor obstacleColor;

		QVector2D cellPixelSize;
		drawModes drawMode;
		double toPixels;
		double uMax;

		Mouse mouse;
		
		Lattice lattice;
};

#endif