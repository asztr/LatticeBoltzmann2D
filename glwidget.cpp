#include "glwidget.h"

const double toDegrees = 180.0 / (4 * atan(1.0));

GLWidget::GLWidget(QWidget *parent) : QGLWidget(QGLFormat(QGL::SampleBuffers), parent) {
	setMinimumSize(640, 480);
	setAutoFillBackground(false);
	grabKeyboard();

	background = QBrush(QColor(0, 85, 255));
	obstacleColor = QColor(0, 0, 0);
	particleColor = QColor(200, 200, 200);
	arrowColor = QColor(255, 255, 0);
	
	arrowPen = QPen(arrowColor);
	particlePen = QPen(particleColor);
	arrowPen.setWidth(1);
	particlePen.setWidth(3);

	double _omega = 0.5;
	double _rho = 1.0;
	int _nParticles = 20000;
	lattice.init(100, 100, false, false, _rho, _omega, _nParticles);
	drawMode = testParticles;
	uMax = lattice.getUMax();
}

void GLWidget::paintParticles(QPainter *painter) {
	painter->setPen(particlePen);
	for (int n=0; n < lattice.getNParticles(); n++) {
		vec<double> _pos = lattice.getParticlePosition(n);
		//painter->drawEllipse(_pos.x*cellPixelSize.x(), _pos.y*cellPixelSize.y(), cellPixelSize.x(), cellPixelSize.y());
		painter->drawPoint(_pos.x*cellPixelSize.x(), _pos.y*cellPixelSize.y());
	}
}

void GLWidget::paintObstacles(QPainter *painter) {
	for (int i = 0; i < lattice.getDim().x; i++)
		for (int j = 0; j < lattice.getDim().y; j++)
			if (lattice.isObstacle(i,j))
				painter->fillRect(floor(i*cellPixelSize.x()), floor(j*cellPixelSize.y()), ceil(cellPixelSize.x()), ceil(cellPixelSize.y()), obstacleColor);
}

void GLWidget::paintVectorField(QPainter *painter) {
	if (drawMode == vectorFieldRenormalized)
		uMax = lattice.getUMax();
	toPixels = cellPixelSize.length()/(sqrt(2.0)*uMax);

	for (int i = 0; i < lattice.getDim().x; i++)
		for (int j = 0; j < lattice.getDim().y; j++)
			if (lattice.isObstacle(i,j) == false) {
				double uMag = lattice.getU(i,j).norm();
				double angle = lattice.getU(i,j).angle() * toDegrees;
				QPointF scale(uMag*toPixels/10, uMag*toPixels/12);
				drawArrow(painter, QPoint(i*cellPixelSize.x(), j*cellPixelSize.y()), scale, angle);
			}
}

void GLWidget::paintFourierField(QPainter *painter) {
	for (int i = 0; i < lattice.getDim().x; i++)
		for (int j = 0; j < lattice.getDim().y; j++)
			if (lattice.isObstacle(i,j) == false) {
				double uMag = lattice.getU(i,j).norm();
				double angle = lattice.getU(i,j).angle() * toDegrees;
				QPointF scale(uMag*toPixels/10, uMag*toPixels/12);
				drawArrow(painter, QPoint(i*cellPixelSize.x(), j*cellPixelSize.y()), scale, angle);
			}
}

void GLWidget::paint(QPainter *painter, QPaintEvent *event) {
	cellPixelSize.setX((double)size().width()/lattice.getDim().x);
	cellPixelSize.setY((double)size().height()/lattice.getDim().y);

	painter->fillRect(event->rect(), background);
	paintObstacles(painter);

	switch(drawMode) {
		case vectorField:
			paintVectorField(painter);
			break;
		case vectorFieldRenormalized:
			paintVectorField(painter);
			break;
		case testParticles:
			paintParticles(painter);
			break;
		case fourierField:
			paintFourierField(painter);
			break;
	}
}

void GLWidget::mouseMoveEvent(QMouseEvent *event) {
	QPoint _pressed;
	if (event->buttons() != Qt::NoButton) {
		mouse.setPosition(event->pos());
		_pressed = mouse.pressedCell(cellPixelSize);
	}
	if (event->buttons() == Qt::LeftButton) {
		QVector2D _delta = mouse.deltaPosition().normalized(); //habria que revisar la eficiencia de crear este objeto
		vec<double> _u(_delta.x(), _delta.y());
		lattice.setU(_pressed.x(), _pressed.y(), _u);
	}
	if (event->buttons() == Qt::RightButton)
		lattice.setObstacle(_pressed.x(), _pressed.y());
}

void GLWidget::keyPressEvent(QKeyEvent *event) {
	if (event->key() == Qt::Key_R) {
		double _rho = 1.0;
		lattice.resetGrid(_rho);
	}
	if (event->key() == Qt::Key_D) {
		drawMode = (drawModes)(1 + ((int)drawMode + 1) % 3);
	}
}

void GLWidget::animate() {
	if (drawMode == testParticles) //este if lo hice a los pedos. quiza no sea la mejor forma de hacerlo, ni el lugar ideal
		lattice.showParticles();
	else
		lattice.hideParticles();
	lattice.update();
	repaint();
}

void GLWidget::paintEvent(QPaintEvent *event) {
	QPainter painter;
	painter.begin(this);
	painter.setRenderHint(QPainter::Antialiasing);
	paint(&painter, event);
	painter.end();
}

void GLWidget::drawArrow(QPainter *painter, QPoint p0, QPointF scale, double angle) {
	painter->setPen(arrowPen);
	painter->save();
	painter->translate(p0);
	painter->rotate(-angle);
	painter->scale(scale.x(), scale.y());
	painter->drawLine(0, 0, 10, 0);
	painter->drawLine(10, 0, 7, -3);
	painter->drawLine(10, 0, 7, 3);
	painter->restore();
}
